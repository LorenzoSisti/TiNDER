import csv

sc_file_path = "/Path/to/the/shape_complementarity/or/electrostatic_complementarity/score_totale.sc"
csv_file_path = sc_file_path.replace(".sc", ".csv")

header_written = False
header = []

with open(sc_file_path, "r") as sc_file, open(csv_file_path, "w", newline="") as csv_file:
    writer = csv.writer(csv_file)

    for line in sc_file:
        if not line.startswith("SCORE:"):
            continue  

        parts = line.strip().split()[1:]  
        
        if 'description' in parts:
            if not header_written:
                header = parts
                writer.writerow(header)
                header_written = True
        else:
            if header_written and len(parts) == len(header):
                writer.writerow(parts)

print(f"✅ Conversion completed: {csv_file_path}")
