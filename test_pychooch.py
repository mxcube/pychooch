import PyChooch
import sys


def main():
    filename = sys.argv[1]
    with open(filename) as f:
        for _ in range(2):
            next(f)
        data = [map(float, line.split()) for line in f]
    result = PyChooch.calc(data, "Se", "K")
    print(result)


if __name__ == '__main__':
    main()
