
import IPDA_analyzer

input_IPDA_file = 'input_files/DCC Titan Plate 5.csv'
input_DNAconcentration = 'input_files/DCC Titan Plate 5_DNAconcentration_input_file.xlsx'
minimum_required_droplets = 10000

run = IPDA_analyzer.export_analyzed_IPDA_as_Excel_and_png(input_IPDA_file, input_DNAconcentration, minimum_required_droplets)
