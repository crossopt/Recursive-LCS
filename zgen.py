n, rep = map(int, input().split())

tries = [""] * (n + 1)
tries[0] = (chr(ord('a') + rep)) * 2
for i in range(1, n + 1):
	tries[i] = tries[i - 1] + (chr(ord('a') + (i + rep) % 26))
fout = open('file', 'w')
print(''.join(tries), file=fout)
fout.close()