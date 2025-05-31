import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--kpoint', type=str, default=None, help='K-point coordinates', required=True)
args = parser.parse_args()


input = open("INCAR", mode='r')

output = open("settings.csv", mode='w')


for line in input:
    parts = line.split("=")

    output.write("{},{}\n".format(parts[0].strip().lower(), parts[1].strip()))

output.write("kpoints,{}\n".format(args.kpoint))