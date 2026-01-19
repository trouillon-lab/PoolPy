import argparse
import os
import re
import time
import numpy as np
import pandas as pd

from Functions_wrapped import *


def parse_method(filename: str) -> str:
	match = re.match(r"^WA_(?P<method>.+)_N_\d+_diff_\d+.*\.csv$", filename)
	return match.group('method') if match else 'unknown'


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--directory', type=str, default='./wrapped', help='Base directory containing N_* folders')
	parser.add_argument('--cleanup', type=str, default='True', help='Remove existing metrics files before computing new ones')
	parser.add_argument('--start', type=int, default=None, help='Starting N value (inclusive)')
	parser.add_argument('--stop', type=int, default=None, help='Stopping N value (inclusive)')
	parser.add_argument('--step', type=int, default=10, help='Step size for N values')
	args = parser.parse_args()
	
	cleanup = args.cleanup in ['True', 'true', True, '1']

	overall_start = time.time()
	total_files = 0
	
	# Determine N values to process
	if args.start is not None and args.stop is not None:
		# Generate N values from start to stop with step
		n_values = range(args.start, args.stop + 1, args.step)
	elif args.start is not None:
		# Start specified but no stop - get all existing folders >= start
		all_n_folders = sorted([f for f in os.listdir(args.directory) if f.startswith('N_')])
		n_values = []
		for n_folder in all_n_folders:
			try:
				n_val = int(n_folder.split('_', 1)[1])
				if n_val >= args.start:
					n_values.append(n_val)
			except ValueError:
				continue
	elif args.stop is not None:
		# Stop specified but no start - get all existing folders <= stop
		all_n_folders = sorted([f for f in os.listdir(args.directory) if f.startswith('N_')])
		n_values = []
		for n_folder in all_n_folders:
			try:
				n_val = int(n_folder.split('_', 1)[1])
				if n_val <= args.stop:
					n_values.append(n_val)
			except ValueError:
				continue
	else:
		# No start/stop specified - get all existing folders
		all_n_folders = sorted([f for f in os.listdir(args.directory) if f.startswith('N_')])
		n_values = []
		for n_folder in all_n_folders:
			try:
				n_values.append(int(n_folder.split('_', 1)[1]))
			except ValueError:
				continue
	
	for n_val in sorted(n_values):
		n_folder = f'N_{n_val}'
		n_dir = os.path.join(args.directory, n_folder)
		if not os.path.isdir(n_dir):
			continue
		
		n_start_time = time.time()
		
		# Get all diff folders for this N
		diff_folders = sorted([f for f in os.listdir(n_dir) if f.startswith('diff_')])
		
		for diff_folder in diff_folders:
			diff_dir = os.path.join(n_dir, diff_folder)
			if not os.path.isdir(diff_dir):
				continue
			
			try:
				diff_val = int(diff_folder.split('_', 1)[1])
			except ValueError:
				continue
			
			diff_start_time = time.time()
			
			# Process WAs folder for this diff
			wa_dir = os.path.join(diff_dir, 'WAs')
			if not os.path.isdir(wa_dir):
				print(f"Skipping {diff_dir}: no WAs folder found")
				continue
			
			# Collect metrics for all WAs in this diff
			rows = []
			wa_files = sorted([f for f in os.listdir(wa_dir) if f.endswith('.csv')])
			
			for wa_file in wa_files:
				wa_path = os.path.join(wa_dir, wa_file)
				method = parse_method(wa_file)
				
				try:
					wa = np.loadtxt(wa_path, delimiter=',', dtype=int)
				except Exception as exc:
					print(f"Skipping {wa_path}: {exc}")
					continue
				
				metrics = compute_metrics(n_val, diff_val, wa)
				metrics['Pooling strategy'] = method
				rows.append(metrics)
				print(f"Processed {method} for N={n_val}, diff={diff_val}")

			# Add Hierarchical method using calculate_metrics_hierarchical_fast with 65 checks
			try:
				hier_metrics = calculate_metrics_hierarchical_fast(n_val, diff_val, checks=65)
				hier_row = {
					'N': n_val,
					'diff': diff_val,
					'Mean experiments': float(np.round(hier_metrics[0], 2)),
					'Max samples per pool': int(hier_metrics[1]),
					'N pools': hier_metrics[2],
					'Percentage check': int(hier_metrics[3]),
					'Max experiments per sample': len(hier_metrics[2])+1,
					'Mean extra experiments': float(np.round(hier_metrics[4], 2)),
					'Mean steps': hier_metrics[5],
					'Pooling strategy': 'Hierarchical',
				}
				rows.append(hier_row)
				print(f"Processed Hierarchical for N={n_val}, diff={diff_val}")
			except Exception as exc:
				print(f"Skipping Hierarchical for N={n_val}, diff={diff_val}: {exc}")
			
			# Save metrics for this N/diff combination
			if rows:
				df = pd.DataFrame(rows)
				numeric_cols = df.select_dtypes(include=[np.number]).columns
				df[numeric_cols] = df[numeric_cols].round(2)
				df.sort_values(['Pooling strategy'], inplace=True)
				
				output_path = os.path.join(diff_dir, f'Metrics_N_{n_val}_diff_{diff_val}.csv')
				
				# Remove existing file if cleanup is enabled
				if cleanup and os.path.exists(output_path):
					os.remove(output_path)
					print(f"Removed existing {output_path}")
				
				df.to_csv(output_path, index=False)
				total_files += 1
				print(f"Saved {output_path} with {len(rows)} entries")
			
			# Print diff completion
			DTS = np.round((time.time() - diff_start_time), 2)
			DTD = DTS // 86400
			DTH = DTS // 3600 - DTD * 24
			DTM = DTS // 60 - DTH * 60 - DTD * 24 * 60
			DTS = np.round(DTS - (DTM + DTH * 60 + DTD * 24 * 60) * 60, 2)
			print("%s days %s hours %s minutes and %s seconds required for N= %s and differentiate %s" % 
					(DTD, DTH, DTM, DTS, n_val, diff_val))
			print('----------------------------------------------------------------------------------------------------------')
		
		# Print N completion
		DTS = np.round((time.time() - n_start_time), 2)
		DTD = DTS // 86400
		DTH = DTS // 3600 - DTD * 24
		DTM = DTS // 60 - DTH * 60 - DTD * 24 * 60
		DTS = np.round(DTS - (DTM + DTH * 60 + DTD * 24 * 60) * 60, 2)
		print('\n')
		print('##########################################################################################################')
		print('##########################################################################################################')
		print('\n')
		print("%s days %s hours %s minutes and %s seconds overall required for N= %s" % (DTD, DTH, DTM, DTS, n_val))
		print('\n')
		print('##########################################################################################################')
		print('##########################################################################################################')
	
	# Print overall completion
	total_time = time.time() - overall_start
	DTS = np.round(total_time, 2)
	DTD = DTS // 86400
	DTH = DTS // 3600 - DTD * 24
	DTM = DTS // 60 - DTH * 60 - DTD * 24 * 60
	DTS = np.round(DTS - (DTM + DTH * 60 + DTD * 24 * 60) * 60, 2)
	print('\n')
	print('\n')
	print('**********************************************************************************************************')
	print('**********************************************************************************************************')
	print('**********************************************************************************************************')
	print('\n')
	print("%s days %s hours %s minutes and %s seconds overall required for run" % (DTD, DTH, DTM, DTS))
	print('\n')
	print(f'Wrote {total_files} CSV files')
	print('\n')
	print('**********************************************************************************************************')
	print('**********************************************************************************************************')
	print('**********************************************************************************************************')


if __name__ == '__main__':
	main()
