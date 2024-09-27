import sys
import statistics


vals = []

for line in sys.stdin:
    num = int(line)
    if num > 0:
        vals.append(num)

print(min(vals))
print(statistics.median(vals))
print(max(vals))
