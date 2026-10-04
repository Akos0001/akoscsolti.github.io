let currentResult = null;

// ============================================================
// Existing calibrations
// ============================================================
//
// Add your already calculated calibrations here.
//
// The values p and cov are the polynomial-fit coefficients and
// covariance matrix returned by calculateEfficiency().
//
// The easiest way to obtain these values is to process the
// original Excel file with this page and copy the corresponding
// values from currentResult.
//
// Example structure:
//
// {
//     id: "hpge_10cm",
//     name: "HPGe - 10 cm",
//     description: "Efficiency calibration at 10 cm source-detector distance.",
//     result: {
//         E: [...],
//         eta: [...],
//         etaAbsUnc: [...],
//         etaRelUnc: [...],
//         p: [...],
//         cov: [[...], [...], [...], [...]],
//         sigmaP: [...],
//         results: [...],
//         Efit: [...],
//         etaFit: [...],
//         sigmaEtaFit: [...],
//         upper: [...],
//         lower: [...]
//     }
// }
//
// IMPORTANT:
// Replace the empty array below with your actual calibrations.
// Once they are added here, visitors can open them without an
// Excel file.

let existingCalibrations = [];


// ============================================================
// Load existing calibrations from JSON files
// ============================================================

async function loadExistingCalibrations() {

    const calibrationFiles = [
        "calibrations/hpge_calibration.json"
    ];

    existingCalibrations = [];

    for (const file of calibrationFiles) {

        try {

            const response = await fetch(
                file + "?v=" + Date.now()
            );

            if (!response.ok) {
                throw new Error(
                    `Could not load ${file}`
                );
            }

            const calibration =
                await response.json();

            existingCalibrations.push(calibration);

        } catch (error) {

            console.error(
                "Error loading calibration:",
                file,
                error
            );
        }
    }
}

let currentExistingCalibration = null;


// ============================================================
// Matrix utilities
// ============================================================

function zeros(rows, cols) {
    return Array.from(
        { length: rows },
        () => Array(cols).fill(0)
    );
}


function identityMatrix(n) {

    const I = zeros(n, n);

    for (let i = 0; i < n; i++) {
        I[i][i] = 1;
    }

    return I;
}


function transpose(A) {

    const rows = A.length;
    const cols = A[0].length;

    const result = zeros(cols, rows);

    for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
            result[j][i] = A[i][j];
        }
    }

    return result;
}


function matrixMultiply(A, B) {

    const rowsA = A.length;
    const colsA = A[0].length;
    const colsB = B[0].length;

    const result = zeros(rowsA, colsB);

    for (let i = 0; i < rowsA; i++) {

        for (let k = 0; k < colsA; k++) {

            const aik = A[i][k];

            for (let j = 0; j < colsB; j++) {
                result[i][j] += aik * B[k][j];
            }
        }
    }

    return result;
}


function matrixVectorMultiply(A, v) {

    return A.map(row =>
        row.reduce(
            (sum, value, i) => sum + value * v[i],
            0
        )
    );
}


function matrixVectorOuterProduct(v) {

    const n = v.length;

    const result = zeros(n, n);

    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            result[i][j] = v[i] * v[j];
        }
    }

    return result;
}


// ============================================================
// Matrix inverse
// ============================================================

function inverseMatrix(A) {

    const n = A.length;

    const M = A.map(
        (row, i) => [
            ...row,
            ...identityMatrix(n)[i]
        ]
    );

    for (let col = 0; col < n; col++) {

        let pivotRow = col;
        let maxValue = Math.abs(M[col][col]);

        for (let row = col + 1; row < n; row++) {

            const value = Math.abs(M[row][col]);

            if (value > maxValue) {
                maxValue = value;
                pivotRow = row;
            }
        }

        if (maxValue < 1e-15) {
            throw new Error("Matrix is singular or nearly singular.");
        }

        if (pivotRow !== col) {
            [M[col], M[pivotRow]] =
                [M[pivotRow], M[col]];
        }

        const pivot = M[col][col];

        for (let j = 0; j < 2 * n; j++) {
            M[col][j] /= pivot;
        }

        for (let row = 0; row < n; row++) {

            if (row === col) {
                continue;
            }

            const factor = M[row][col];

            for (let j = 0; j < 2 * n; j++) {
                M[row][j] -= factor * M[col][j];
            }
        }
    }

    return M.map(
        row => row.slice(n)
    );
}


// ============================================================
// QR decomposition using modified Gram-Schmidt
// ============================================================

function qrDecomposition(A) {

    const m = A.length;
    const n = A[0].length;

    const Q = zeros(m, n);
    const R = zeros(n, n);

    for (let j = 0; j < n; j++) {

        let v = A.map(row => row[j]);

        for (let i = 0; i < j; i++) {

            let dot = 0;

            for (let k = 0; k < m; k++) {
                dot += Q[k][i] * v[k];
            }

            R[i][j] = dot;

            for (let k = 0; k < m; k++) {
                v[k] -= dot * Q[k][i];
            }
        }

        let norm = Math.sqrt(
            v.reduce(
                (sum, value) => sum + value * value,
                0
            )
        );

        if (norm < 1e-14) {
            throw new Error(
                "The weighted fitting matrix is rank deficient."
            );
        }

        R[j][j] = norm;

        for (let k = 0; k < m; k++) {
            Q[k][j] = v[k] / norm;
        }
    }

    return { Q, R };
}


// ============================================================
// Solve upper triangular system R*x = b
// ============================================================

function solveUpperTriangular(R, b) {

    const n = R.length;
    const x = new Array(n).fill(0);

    for (let i = n - 1; i >= 0; i--) {

        let sum = b[i];

        for (let j = i + 1; j < n; j++) {
            sum -= R[i][j] * x[j];
        }

        x[i] = sum / R[i][i];
    }

    return x;
}


// ============================================================
// Weighted polynomial fitting
// ============================================================

function weightedPolynomialFit(logE, logEta, sigmaLogEta) {

    const n = logE.length;
    const deg = 3;

    const A = zeros(n, deg + 1);

    for (let i = 0; i < n; i++) {

        const x = logE[i];

        A[i][0] = x ** 3;
        A[i][1] = x ** 2;
        A[i][2] = x;
        A[i][3] = 1;
    }

    const Aw = zeros(n, deg + 1);
    const bw = new Array(n);

    for (let i = 0; i < n; i++) {

        const weight = 1 / sigmaLogEta[i];

        for (let j = 0; j < deg + 1; j++) {
            Aw[i][j] = weight * A[i][j];
        }

        bw[i] = weight * logEta[i];
    }

    const { Q, R } = qrDecomposition(Aw);

    const Qt = transpose(Q);

    const Qtb = matrixVectorMultiply(Qt, bw);

    const p = solveUpperTriangular(R, Qtb);

    const Rinv = inverseMatrix(R);

    const cov = matrixMultiply(
        Rinv,
        transpose(Rinv)
    );

    const sigmaP = [];

    for (let i = 0; i < cov.length; i++) {

        sigmaP.push(
            Math.sqrt(
                Math.max(cov[i][i], 0)
            )
        );
    }

    return {
        p,
        cov,
        sigmaP
    };
}


// ============================================================
// Polynomial evaluation
// ============================================================

function polynomialValue(p, x) {

    return (
        p[0] * x ** 3 +
        p[1] * x ** 2 +
        p[2] * x +
        p[3]
    );
}


// ============================================================
// Efficiency calculation
// ============================================================

function calculateEfficiency(rows) {

    if (rows.length === 0) {
        throw new Error("The Excel file contains no data.");
    }

    const E = [];
    const k = [];
    const N = [];
    const NRelUnc = [];
    const tLive = [];
    const tRef = [];
    const A0 = [];
    const A0Unc = [];
    const T = [];

    for (const row of rows) {

        E.push(Number(row[0]));

        k.push(Number(row[1]) * 1e-2);

        N.push(Number(row[2]));

        NRelUnc.push(Number(row[3]) * 1e-2);

        tLive.push(Number(row[4]));

        tRef.push(Number(row[7]));

        A0.push(Number(row[8]) * 1e3);

        A0Unc.push(Number(row[9]) * 1e3);

        T.push(Number(row[10]));
    }

    for (let i = 0; i < E.length; i++) {

        const values = [
            E[i],
            k[i],
            N[i],
            NRelUnc[i],
            tLive[i],
            tRef[i],
            A0[i],
            A0Unc[i],
            T[i]
        ];

        if (values.some(v => !Number.isFinite(v))) {

            throw new Error(
                `Invalid numerical data in row ${i + 2}.`
            );
        }
    }

    const A0RelUnc = [];
    const A = [];
    const ARelUnc = [];
    const AAbsUnc = [];

    for (let i = 0; i < E.length; i++) {

        const relUnc = A0Unc[i] / A0[i];

        const activity =
            A0[i] *
            2 ** (-tRef[i] / T[i]);

        const absUnc =
            relUnc * activity;

        A0RelUnc.push(relUnc);
        A.push(activity);
        ARelUnc.push(relUnc);
        AAbsUnc.push(absUnc);
    }

    const eta = [];
    const etaRelUnc = [];
    const etaAbsUnc = [];

    for (let i = 0; i < E.length; i++) {

        const efficiency =
            N[i] /
            (A[i] * tLive[i] * k[i]);

        const relativeUnc =
            Math.sqrt(
                NRelUnc[i] ** 2 +
                ARelUnc[i] ** 2
            );

        const absoluteUnc =
            efficiency * relativeUnc;

        eta.push(efficiency);
        etaRelUnc.push(relativeUnc);
        etaAbsUnc.push(absoluteUnc);
    }

    const EeV = E.map(
        value => value * 1e3
    );

    const logE = EeV.map(
        value => Math.log(value)
    );

    const logEta = eta.map(
        value => Math.log(value)
    );

    const sigmaLogEta = etaAbsUnc.map(
        (unc, i) => unc / eta[i]
    );

    const fit = weightedPolynomialFit(
        logE,
        logEta,
        sigmaLogEta
    );

    const p = fit.p;
    const cov = fit.cov;
    const sigmaP = fit.sigmaP;

    const etaFitData = [];

    for (let i = 0; i < logE.length; i++) {

        const logEtaFit =
            polynomialValue(
                p,
                logE[i]
            );

        etaFitData.push(
            Math.exp(logEtaFit)
        );
    }

    const absDiff = [];
    const relDiff = [];

    for (let i = 0; i < eta.length; i++) {

        const difference =
            eta[i] - etaFitData[i];

        absDiff.push(difference);

        relDiff.push(
            100 * difference / etaFitData[i]
        );
    }

    const minE =
        Math.min(...EeV);

    const maxE =
        Math.max(...EeV);

    const Efit = [];

    const nFit = 1000;

    const logMin = Math.log(0.8 * minE);
    const logMax = Math.log(1.2 * maxE);

    for (let i = 0; i < nFit; i++) {

        const x =
            logMin +
            i * (logMax - logMin) / (nFit - 1);

        Efit.push(
            Math.exp(x)
        );
    }

    const etaFit = [];
    const sigmaEtaFit = [];
    const upper = [];
    const lower = [];

    for (const energy of Efit) {

        const x = Math.log(energy);

        const logEtaFit =
            polynomialValue(p, x);

        const efficiency =
            Math.exp(logEtaFit);

        const J = [
            x ** 3,
            x ** 2,
            x,
            1
        ];

        const covJ =
            matrixVectorMultiply(
                cov,
                J
            );

        let variance = 0;

        for (let i = 0; i < J.length; i++) {
            variance +=
                J[i] * covJ[i];
        }

        const sigmaLogEta =
            Math.sqrt(
                Math.max(variance, 0)
            );

        const sigmaEta =
            efficiency * sigmaLogEta;

        etaFit.push(efficiency);
        sigmaEtaFit.push(sigmaEta);

        upper.push(
            efficiency + sigmaEta
        );

        lower.push(
            Math.max(
                efficiency - sigmaEta,
                1e-10
            )
        );
    }

    const results = [];

    for (let i = 0; i < E.length; i++) {

        results.push({
            energyKeV: EeV[i] * 1e-3,
            measuredEfficiency: eta[i],
            fittedEfficiency: etaFitData[i],
            absoluteDifference: absDiff[i],
            relativeDifference: relDiff[i]
        });
    }

    return {
        E: EeV,
        eta,
        etaAbsUnc,
        etaRelUnc,

        p,
        cov,
        sigmaP,

        results,

        Efit,
        etaFit,
        sigmaEtaFit,
        upper,
        lower
    };
}


// ============================================================
// Efficiency at given energy
// ============================================================

function calculateEfficiencyAtEnergy(
    EkeV,
    p,
    cov
) {

    const EeV =
        EkeV * 1e3;

    const x =
        Math.log(EeV);

    const J = [
        x ** 3,
        x ** 2,
        x,
        1
    ];

    const logEta =
        J.reduce(
            (sum, value, i) =>
                sum + value * p[i],
            0
        );

    const eta =
        Math.exp(logEta);

    const covJ =
        matrixVectorMultiply(
            cov,
            J
        );

    let variance = 0;

    for (let i = 0; i < J.length; i++) {

        variance +=
            J[i] * covJ[i];
    }

    const sigmaLogEta =
        Math.sqrt(
            Math.max(variance, 0)
        );

    const sigmaEta =
        eta * sigmaLogEta;

    const relativeUncertainty =
        100 * sigmaEta / eta;

    return {
        energyKeV: EkeV,
        eta,
        sigmaEta,
        relativeUncertainty
    };
}


// ============================================================
// Plot
// ============================================================

function createPlot(result, plotId = "plot") {

    const measuredTrace = {
        x: result.E,
        y: result.eta,
        error_y: {
            type: "data",
            array: result.etaAbsUnc,
            visible: true
        },
        mode: "markers",
        name: "Measured data",
        marker: {
            size: 7
        }
    };

    const fitTrace = {
        x: result.Efit,
        y: result.etaFit,
        mode: "lines",
        name: "Fitted curve",
        line: {
            width: 2
        }
    };

    const upperTrace = {
        x: result.Efit,
        y: result.upper,
        mode: "lines",
        name: "Fit uncertainty",
        line: {
            dash: "dash",
            width: 1
        },
        showlegend: false
    };

    const lowerTrace = {
        x: result.Efit,
        y: result.lower,
        mode: "lines",
        line: {
            dash: "dash",
            width: 1
        },
        showlegend: false
    };

    const uncertaintyBand = {
        x: [
            ...result.Efit,
            ...result.Efit.slice().reverse()
        ],

        y: [
            ...result.upper,
            ...result.lower.slice().reverse()
        ],

        fill: "toself",
        fillcolor: "rgba(100,100,100,0.15)",
        line: {
            color: "rgba(255,255,255,0)"
        },
        name: "Fit uncertainty",
        hoverinfo: "skip",
        showlegend: false
    };

    const layout = {

        title:
            "ln(ε) = a + b ln(E) + c (ln(E))² + d (ln(E))³",

        xaxis: {
            title: "Energy [eV]",
            type: "log"
        },

        yaxis: {
            title: "Efficiency [-]",
            type: "log"
        },

        hovermode: "closest",

        legend: {
            orientation: "h"
        },

        margin: {
            l: 70,
            r: 30,
            t: 80,
            b: 70
        }
    };

    Plotly.newPlot(
        plotId,
        [
            uncertaintyBand,
            fitTrace,
            upperTrace,
            lowerTrace,
            measuredTrace
        ],
        layout,
        {
            responsive: true
        }
    );
}


// ============================================================
// Results table
// ============================================================

function updateResultsTable(results) {

    const tbody =
        document.querySelector(
            "#resultsTable tbody"
        );

    tbody.innerHTML = "";

    for (const row of results) {

        const tr =
            document.createElement("tr");

        tr.innerHTML = `
            <td>${row.energyKeV.toPrecision(6)}</td>
            <td>${row.measuredEfficiency.toExponential(6)}</td>
            <td>${row.fittedEfficiency.toExponential(6)}</td>
            <td>${row.absoluteDifference.toExponential(6)}</td>
            <td>${row.relativeDifference.toFixed(4)}</td>
        `;

        tbody.appendChild(tr);
    }
}


// ============================================================
// Message
// ============================================================

function showMessage(
    text,
    type
) {

    const message =
        document.getElementById("message");

    message.textContent = text;

    message.className =
        `message ${type}`;
}


// ============================================================
// Excel processing
// ============================================================

async function processExcelFile(file) {

    if (!file) {

        showMessage(
            "Please choose an Excel file first.",
            "error"
        );

        return;
    }

    try {

        showMessage(
            "Processing file...",
            "success"
        );

        const arrayBuffer =
            await file.arrayBuffer();

        const workbook =
            XLSX.read(
                arrayBuffer,
                {
                    type: "array"
                }
            );

        const firstSheetName =
            workbook.SheetNames[0];

        const worksheet =
            workbook.Sheets[firstSheetName];

        const rows =
            XLSX.utils.sheet_to_json(
                worksheet,
                {
                    header: 1,
                    defval: null
                }
            );

        const dataRows =
            rows.slice(1).filter(
                row =>
                    row.length > 0 &&
                    row[0] !== null &&
                    row[0] !== ""
            );

        currentResult =
            calculateEfficiency(dataRows);
        console.log("CALIBRATION RESULT:");
        console.log(JSON.stringify(currentResult, null, 4));

        createPlot(
            currentResult
        );

        updateResultsTable(
            currentResult.results
        );

        document.getElementById(
            "resultsCard"
        ).style.display = "block";

        showMessage(
            `Successful processing: ${file.name}`,
            "success"
        );

    } catch (error) {

        console.error(error);

        showMessage(
            "Error during processing: " +
            error.message,
            "error"
        );

        document.getElementById(
            "resultsCard"
        ).style.display = "none";
    }
}


// ============================================================
// Existing calibration interface
// ============================================================

function showPage(page) {

    document.getElementById(
        "startPage"
    ).classList.add("section-hidden");

    document.getElementById(
        "newCalibrationPage"
    ).classList.add("section-hidden");

    document.getElementById(
        "existingCalibrationPage"
    ).classList.add("section-hidden");

    if (page === "start") {

        document.getElementById(
            "startPage"
        ).classList.remove("section-hidden");
    }

    if (page === "new") {

        document.getElementById(
            "newCalibrationPage"
        ).classList.remove("section-hidden");
    }

    if (page === "existing") {

        document.getElementById(
            "existingCalibrationPage"
        ).classList.remove("section-hidden");
    }
}


function showExistingMessage(
    text,
    type
) {

    const message =
        document.getElementById(
            "existingCalibrationMessage"
        );

    message.textContent = text;

    message.className =
        `message ${type}`;
}


function renderExistingCalibrationList() {

    const list =
        document.getElementById(
            "calibrationList"
        );

    list.innerHTML = "";

    if (existingCalibrations.length === 0) {

        const empty =
            document.createElement("div");

        empty.className =
            "calibration-card";

        empty.innerHTML = `
            <div>
                <h3>No existing calibrations available</h3>
                <p>
                    Add your calculated calibrations to the
                    <code>existingCalibrations</code> array in script.js.
                </p>
            </div>
        `;

        list.appendChild(empty);

        return;
    }

    for (const calibration of existingCalibrations) {

        const card =
            document.createElement("div");

        card.className =
            "calibration-card";

        card.innerHTML = `
            <div>
                <h3>${calibration.name}</h3>
                <p>${calibration.description || ""}</p>
            </div>

            <button
                class="button"
                type="button"
                data-calibration-id="${calibration.id}"
            >
                Open
            </button>
        `;

        card.querySelector("button").addEventListener(
            "click",
            () => openExistingCalibration(
                calibration.id
            )
        );

        list.appendChild(card);
    }
}


function updateExistingResultsTable(results) {

    const tbody =
        document.querySelector(
            "#existingResultsTable tbody"
        );

    tbody.innerHTML = "";

    for (const row of results) {

        const tr =
            document.createElement("tr");

        tr.innerHTML = `
            <td>${row.energyKeV.toPrecision(6)}</td>
            <td>${row.measuredEfficiency.toExponential(6)}</td>
            <td>${row.fittedEfficiency.toExponential(6)}</td>
            <td>${row.absoluteDifference.toExponential(6)}</td>
            <td>${row.relativeDifference.toFixed(4)}</td>
        `;

        tbody.appendChild(tr);
    }
}


function openExistingCalibration(id) {

    const calibration =
        existingCalibrations.find(
            item => item.id === id
        );

    if (
        !calibration ||
        !calibration.result
    ) {

        showExistingMessage(
            "The selected calibration is not available or is incomplete.",
            "error"
        );

        return;
    }

    currentExistingCalibration =
        calibration;

    const result =
        calibration.result;

    document.getElementById(
        "existingCalibrationTitle"
    ).textContent =
        calibration.name;

    document.getElementById(
        "existingCalibrationInfo"
    ).innerHTML = `
        <p>${calibration.description || ""}</p>
    `;

    createPlot(
        result,
        "existingPlot"
    );

    updateExistingResultsTable(
        result.results
    );

    document.getElementById(
        "existingResultsCard"
    ).style.display =
        "block";

    document.getElementById(
        "existingEnergyResult"
    ).style.display =
        "none";

    showExistingMessage(
        `Calibration opened: ${calibration.name}`,
        "success"
    );

    document.getElementById(
        "existingResultsCard"
    ).scrollIntoView({
        behavior: "smooth",
        block: "start"
    });
}


function calculateExistingEfficiency() {

    if (!currentExistingCalibration) {

        showExistingMessage(
            "Please open a calibration first.",
            "error"
        );

        return;
    }

    const energy =
        Number(
            document.getElementById(
                "existingEnergyInput"
            ).value
        );

    if (
        !Number.isFinite(energy) ||
        energy <= 0
    ) {

        showExistingMessage(
            "Energy must be a positive number.",
            "error"
        );

        return;
    }

    try {

        const result =
            calculateEfficiencyAtEnergy(
                energy,
                currentExistingCalibration.result.p,
                currentExistingCalibration.result.cov
            );

        document.getElementById(
            "existingResultEnergy"
        ).textContent =
            result.energyKeV.toFixed(6) +
            " keV";

        document.getElementById(
            "existingResultEfficiency"
        ).textContent =
            result.eta.toExponential(6);

        document.getElementById(
            "existingResultUncertainty"
        ).textContent =
            "± " +
            result.sigmaEta.toExponential(6);

        document.getElementById(
            "existingResultRelativeUncertainty"
        ).textContent =
            result.relativeUncertainty.toFixed(4) +
            " %";

        document.getElementById(
            "existingEnergyResult"
        ).style.display =
            "block";

    } catch (error) {

        showExistingMessage(
            "Error during efficiency calculation: " +
            error.message,
            "error"
        );
    }
}


// ============================================================
// Event listeners
// ============================================================

document.addEventListener(
    "DOMContentLoaded",
    async () => {

        await loadExistingCalibrations();

        // --------------------------------------------------------
        // Start menu
        // --------------------------------------------------------

        document.getElementById(
            "newCalibrationButton"
        ).addEventListener(
            "click",
            () => {
                showPage("new");
            }
        );

        document.getElementById(
            "existingCalibrationButton"
        ).addEventListener(
            "click",
            () => {

                showPage("existing");

                renderExistingCalibrationList();
            }
        );

        document.getElementById(
            "backToStartFromNew"
        ).addEventListener(
            "click",
            () => {
                showPage("start");
            }
        );

        document.getElementById(
            "backToStartFromExisting"
        ).addEventListener(
            "click",
            () => {
                showPage("start");
            }
        );

        document.getElementById(
            "existingCalculateButton"
        ).addEventListener(
            "click",
            calculateExistingEfficiency
        );

        document.getElementById(
            "existingEnergyInput"
        ).addEventListener(
            "keydown",
            event => {

                if (event.key === "Enter") {
                    calculateExistingEfficiency();
                }
            }
        );

        const fileInput =
            document.getElementById(
                "fileInput"
            );

        const fileName =
            document.getElementById(
                "fileName"
            );

        const processButton =
            document.getElementById(
                "processButton"
            );

        const calculateButton =
            document.getElementById(
                "calculateButton"
            );
        const saveCalibrationButton =
            document.getElementById(
                "saveCalibrationButton"
            );

        const energyInput =
            document.getElementById(
                "energyInput"
            );


        // --------------------------------------------------------
        // File selection
        // --------------------------------------------------------

        fileInput.addEventListener(
            "change",
            () => {

                if (
                    fileInput.files &&
                    fileInput.files.length > 0
                ) {

                    fileName.textContent =
                        fileInput.files[0].name;

                } else {

                    fileName.textContent =
                        "No file selected";
                }
            }
        );


        // --------------------------------------------------------
        // Process button
        // --------------------------------------------------------

        processButton.addEventListener(
            "click",
            () => {

                const file =
                    fileInput.files[0];

                processExcelFile(file);
            }
        );

        // --------------------------------------------------------
        // Save calibration
        // --------------------------------------------------------
        
        saveCalibrationButton.addEventListener(
            "click",
            () => {
        
                if (!currentResult) {
        
                    showMessage(
                        "Please process an Excel file first.",
                        "error"
                    );
        
                    return;
                }
        
                const calibrationData = {
                    id: "hpge_calibration",
                    name: "HPGe detector calibration",
                    description: "Efficiency calibration",
                    result: currentResult
                };
        
                const json =
                    JSON.stringify(
                        calibrationData,
                        null,
                        4
                    );
        
                const blob =
                    new Blob(
                        [json],
                        {
                            type: "application/json"
                        }
                    );
        
                const url =
                    URL.createObjectURL(blob);
        
                const link =
                    document.createElement("a");
        
                link.href = url;
                link.download =
                    "hpge_calibration.json";
        
                document.body.appendChild(link);
        
                link.click();
        
                document.body.removeChild(link);
        
                URL.revokeObjectURL(url);
        
                showMessage(
                    "Calibration saved successfully.",
                    "success"
                );
            }
        );


        // --------------------------------------------------------
        // Energy calculation
        // --------------------------------------------------------

        calculateButton.addEventListener(
            "click",
            () => {

                if (!currentResult) {

                    showMessage(
                        "Please process an Excel file first.",
                        "error"
                    );

                    return;
                }

                const energy =
                    Number(
                        energyInput.value
                    );

                if (
                    !Number.isFinite(energy) ||
                    energy <= 0
                ) {

                    showMessage(
                        "Energy must be a positive number.",
                        "error"
                    );

                    return;
                }

                try {

                    const result =
                        calculateEfficiencyAtEnergy(
                            energy,
                            currentResult.p,
                            currentResult.cov
                        );

                    document.getElementById(
                        "resultEnergy"
                    ).textContent =
                        result.energyKeV.toFixed(6) +
                        " keV";

                    document.getElementById(
                        "resultEfficiency"
                    ).textContent =
                        result.eta.toExponential(6);

                    document.getElementById(
                        "resultUncertainty"
                    ).textContent =
                        "± " +
                        result.sigmaEta.toExponential(6);

                    document.getElementById(
                        "resultRelativeUncertainty"
                    ).textContent =
                        result.relativeUncertainty.toFixed(4) +
                        " %";

                    document.getElementById(
                        "energyResult"
                    ).style.display =
                        "block";

                } catch (error) {

                    showMessage(
                        "Error during efficiency calculation: " +
                        error.message,
                        "error"
                    );
                }
            }
        );


        // --------------------------------------------------------
        // Allow Enter in energy field
        // --------------------------------------------------------

        energyInput.addEventListener(
            "keydown",
            event => {

                if (event.key === "Enter") {
                    calculateButton.click();
                }
            }
        );
    }
);

