from flask import Flask, render_template, request
from calculation import calculate_efficiency, calculate_efficiency_at_energy
import os

app = Flask(__name__)

UPLOAD_FOLDER = "uploads"
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER

os.makedirs(UPLOAD_FOLDER, exist_ok=True)

current_result = None


@app.route("/", methods=["GET", "POST"])
def home():
    global current_result

    message = ""
    energy_result = None

    # Excel feldolgozása
    file = request.files.get("file")

    if file and file.filename:
        filepath = os.path.join(app.config["UPLOAD_FOLDER"], file.filename)
        file.save(filepath)

        try:
            current_result = calculate_efficiency(filepath)
            message = f"Successful processing: {file.filename}"
        except Exception as e:
            current_result = None
            message = f"Hiba a feldolgozás során: {e}"

    # Hatásfok számítása adott energián
    energy = request.form.get("energy")

    if energy:
        try:
            E_keV = float(energy)

            if current_result is not None:
                energy_result = calculate_efficiency_at_energy(
                    E_keV,
                    current_result["p"],
                    current_result["cov"]
                )
            else:
                message = "Először tölts fel és dolgozz fel egy Excel fájlt."

        except ValueError:
            message = "Az energiát számként kell megadni."
        except Exception as e:
            message = f"Hiba a hatásfok számítása során: {e}"

    return render_template(
        "index.html",
        message=message,
        result=current_result,
        energy_result=energy_result
    )


if __name__ == "__main__":
    app.run(debug=True)
