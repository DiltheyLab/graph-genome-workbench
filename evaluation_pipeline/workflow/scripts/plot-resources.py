import sys
from matplotlib.backends.backend_pdf import PdfPages
import gzip
import matplotlib.pyplot as plt
import argparse


def extract_resources(filename):
	cpu_time = 0.0
	highest_max_rss = 0.0
	wallclock_time = 0.0
	for line in open(filename, 'r'):
		fields = line.strip().split()
		if line.strip().startswith("User time (seconds): "):
			cpu_time += float(fields[-1])
		elif line.strip().startswith("System time (seconds):"):
			cpu_time += float(fields[-1])
		elif line.strip().startswith("Maximum resident set size (kbytes):"):
			max_rss = int(fields[-1])
			if highest_max_rss < max_rss:
				highest_max_rss = max_rss
		elif line.strip().startswith("Elapsed (wall clock) time"):
			times = [f for f in fields[-1].split(':')]
			wallclock_time += int(float(times[-1])) + 60.0 * int(float(times[-2]))
			if len(times) > 2:
				wallclock_time += int(float(times[-3])) * 3600.0
		else:
			continue
	return wallclock_time, highest_max_rss, cpu_time


def plot_resources(files, outname, samples, versions):
	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080', '#f7ce37', '#bc202c', '#251188' ]
	runtimes = {}
	wallclock_times = {}
	rss = {}

	for file in files:
		sample = file.strip().split('/')[-4]
		version = file.strip().split('/')[-3]
		wallclock_time, max_rss, cpu_time = extract_resources(file)
		runtimes[(sample, version)] = cpu_time
		rss[(sample, version)] = max_rss
		wallclock_times[(sample, version)] = wallclock_time

	x_values = [i for i in range(len(samples))]
	
	with PdfPages(outname + '.pdf') as pdf:
		# plot runtimes
		fig, ax = plt.subplots()
		for version in versions:
			line_runtimes = [runtimes[(sample, version)] for sample in samples]
			ax.plot(x_values, line_runtimes, label=version, marker='o')
		ax.set_xticks(x_values)
		ax.set_xticklabels(samples, rotation='vertical')
		ax.set_ylabel('Single core CPU seconds')
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# plot wallclock times
		fig, ax = plt.subplots()
		for version in versions:
			line_wallclock = [wallclock_times[(sample, version)] for sample in samples]
			ax.plot(x_values, line_wallclock, label=version, marker='o')
		ax.set_xticks(x_values)
		ax.set_xticklabels(samples, rotation='vertical')
		ax.set_ylabel('Wallclock seconds')
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		plt.tight_layout()
		pdf.savefig()
		plt.close()	

		# plot resources
		fig, ax = plt.subplots()
		for version in versions:
			line_rss = [ rss[(sample, version)]/10**6 for sample in samples]
			ax.plot(x_values, line_rss, label=version, marker='o')
		ax.set_xticks(x_values)
		ax.set_xticklabels(samples, rotation='vertical')
		ax.set_ylabel('Max RSS [GB]')
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		plt.tight_layout()
		pdf.savefig()
		plt.close()



parser = argparse.ArgumentParser(prog='plot-resources.py', description="Plot resources.")
parser.add_argument('-files', metavar='FILES', nargs='+', help='files with results per sample.')
parser.add_argument('-outname', metavar='OUTNAME', required=True, help='Name of the output file.')
parser.add_argument('-samples', metavar='SAMPLES', nargs='+', required=True, help='samples to consider.')
parser.add_argument('-versions', metavar='VERSIONS', nargs='+', required=True, help='Algorithms for genotyping.')
args = parser.parse_args()

plot_resources(args.files, args.outname, [s for s in args.samples], [s for s in args.versions])
