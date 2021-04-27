#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def plot_graph(xs, x2s, ys, method_type):
	plt.plot(ys, x2s, label='Наш алгоритм')
	plt.plot(ys, xs, label='Динамика')
	plt.title('{}'.format(method_type))
	plt.xlabel('Размер сжатого текста (Кб)')
	plt.ylabel('Среднее время работы (сек)')
	plt.legend()
	plt.show()


def plot_agrep_graph(xs, x1s, x2s, ys, method_type):
	plt.plot(ys, x1s, label='Наш алгоритм')
	plt.plot(ys, xs, label='Динамика')
	plt.plot(ys, x2s, label='agrep', color='red')	
	plt.title('{}'.format(method_type))
	plt.xlabel('Размер сжатого текста (Кб)')
	plt.ylabel('Среднее время работы (сек)')
	plt.legend()
	plt.show()



def read_results(file_name):
	with open(file_name, 'r') as fin:
		dps, lzs = [], []
		sizes = []
		for line in fin.readlines():
			p, s, l, d = line.split('&')
			dps.append(float(d) / 1000)
			lzs.append(float(l) / 1000)
			sizes.append(float(s) / 1000)
		return sizes, dps, lzs


if __name__ == "__main__":
	sizes, dps, lzs = read_results('graph3')
	plot_graph(dps, lzs, sizes, "Время подсчета НОП для LZW-сжатия")
	# plot_agrep_graph(dps, lzs, [i * 2 for i in lzs], sizes, "Время подсчета НОП для сжатия UNIX-compress")
