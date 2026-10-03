from calculation import calculate_efficiency


filename = r"C:\Users\akosc\Desktop\Uni\Efficiency_cal_webapp\uploads\meresi_eredmenyek_12356kupak.xlsx"
result = calculate_efficiency(filename)

print()
print("Efficiency results:")
print()

print(result["results"])

print()
print("Fitting coefficients:")
print(result["p"])

print()
print("Coefficient uncertainties:")
print(result["sigma_p"])