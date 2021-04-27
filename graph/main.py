#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def plot_graph(xs, x2s, ys, method_type):
	plt.plot(ys, xs, label='Динамика')
	plt.plot(ys, x2s, label='Наш алгоритм')
	plt.title('{}'.format(method_type))
	plt.xlabel('Размер сжатого текста (байты)')
	plt.ylabel('Среднее время работы (мс)')
	plt.legend()
	plt.show()


def read_results(file_name):
	with open(file_name, 'r') as fin:
		dps, lzs = [], []
		sizes = []
		for line in fin.readlines():
			p, s, l, d = line.split('&')
			dps.append(float(d))
			lzs.append(float(l))
			sizes.append(int(s))
		return sizes, dps, lzs


if __name__ == "__main__":
	sizes, dps, lzs = read_results('graph2')
	plot_graph(dps, lzs, sizes, "Время подсчета НОП для LZW-сжатия")
