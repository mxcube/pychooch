import PyChooch
import sys

if __name__=='__main__':
  filename = sys.argv[1]
  f = open(filename)
  data=[]
  i=0

  for l in f:
    i=i+1
    if (i>2):
      data.append(map(float, l.split()))

  print data

  print PyChooch.calc(data, "Se", "K")
   
