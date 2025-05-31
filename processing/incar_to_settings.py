input = open("INCAR", mode='r')

output = open("settings.csv", mode='w')


for line in input:
    parts = line.split("=")

    output.write("{},{}\n".format(parts[0].strip(), parts[1].strip()))
