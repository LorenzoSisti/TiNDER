import csv

sc_file_path = "/Users/lorenzosisti/Downloads/interface_descriptors/shape_complementarity/score_totale.sc"
csv_file_path = sc_file_path.replace(".sc", ".csv")

header_written = False
header = []

with open(sc_file_path, "r") as sc_file, open(csv_file_path, "w", newline="") as csv_file:
    writer = csv.writer(csv_file)

    for line in sc_file:
        if not line.startswith("SCORE:"):
            continue  # ignora righe che non iniziano con SCORE:

        parts = line.strip().split()[1:]  # rimuove "SCORE:"
        
        # Se è una riga di header (contiene 'description'), salva come intestazione
        if 'description' in parts:
            if not header_written:
                header = parts
                writer.writerow(header)
                header_written = True
        else:
            # Se è una riga di dati e l'header è già stato scritto, scrivi la riga
            if header_written and len(parts) == len(header):
                writer.writerow(parts)

print(f"✅ Conversione completata: {csv_file_path}")
