#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def plot_graph(xs, x2s, ys, method_type):
	fs = 17
	plt.plot(ys, x2s, label='Наш алгоритм')
	plt.plot(ys, xs, label='Динамика')
	plt.title('{}'.format(method_type), fontsize=fs + 1)
	plt.xlabel('Размер сжатого текста (Кб)', fontsize=fs)
	plt.ylabel('Среднее время работы (сек)', fontsize=fs)
	plt.legend(prop={'size': fs})
	plt.show()


def plot_agrep_graph(xs, x1s, x2s, ys, method_type):
	plt.figure(figsize=(8,4))	
	plt.plot(ys, x1s, label='Наш алгоритм')
	plt.plot(ys, xs, label='Динамика')
	plt.plot(ys, x2s, label='agrep', color='red')	
	plt.title('{}'.format(method_type))
	plt.xlabel('Размер сжатого текста (Кб)', fontsize=14)
	plt.ylabel('Среднее время работы (сек)', fontsize=14)
	plt.legend(prop={'size': 14})
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


def read_agrep_results(file_name):
	with open(file_name, 'r') as fin:
		dps, lzs, agrs = [], [], []
		sizes = []
		for line in fin.readlines():
			x, p, s, l, d = line.split(' ')
			dps.append(int(d) // 200000000)
			lzs.append(int(l) // 200000000)
			agrs.append(int(s) // 200000000)
			sizes.append(4 * float(p) / 1000)
		return sizes, dps, lzs, agrs


if __name__ == "__main__":
	#sizes, dps, lzs, agrs = read_agrep_results('graph6')
	sizes, dps, lzs = read_results('lzwreal')
	plot_graph(dps, lzs, sizes, "Время подсчета НОП для LZW-сжатия на естественных текстах")
	# plot_agrep_graph(dps, lzs, agrs, sizes, "Время подсчета НОП для сжатия UNIX-compress")
