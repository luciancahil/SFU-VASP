import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--kpoint', type=str, default=None, help='K-point coordinates', required=True)
args = parser.parse_args()


input = open("Ru2Bi2O7.txt", mode='r')

output = open("Ru2Bi2O7.csv", mode='x')


for line in input:
    parts = line.split("=")

    if(parts[0].strip() == "MAGMOM"):
        output.write("{},{}\n".format(parts[0].strip().lower(), -1))
        print("WARNING: Make sure to edit your MAGMOM field according to the README. The file will not work currently")
    else:
        output.write("{},{}\n".format(parts[0].strip().lower(), parts[1].strip()))

output.write("kpts,{}\n".format(args.kpoint))


