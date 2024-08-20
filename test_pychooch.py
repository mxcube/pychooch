import PyChooch
import sys


def main():
    filename = sys.argv[1]
    with open(filename) as f:
        # Skip first 2 lines
        f.readline()
        f.readline()
        data = [tuple(map(float, line.split())) for line in f]
    result = PyChooch.calc(data, "Se", "K")
    print(result)


if __name__ == '__main__':
    main()
