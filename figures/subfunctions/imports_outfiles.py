#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:18:29 2024

@author: marine
"""
import numpy as np
import sys

class FileReader:
    """A class for reading and processing various output files."""

    @staticmethod
    def read_single_column_file(filename):
        """Reads a file with a single column of data."""
        try:
            with open(filename, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            print(f"Error: {filename} does not exist!")
            sys.exit()
        return np.array([float(line.split()[0]) for line in lines])

    @staticmethod
    def read_dataset(filename):
        """Reads a dataset with time and magnitude columns."""
        try:
            with open(filename, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            print(f"Error: {filename} does not exist!")
            sys.exit()
        times, magnitudes = zip(*[(float(line.split()[0]), float(line.split()[1])) for line in lines])
        return list(times), list(magnitudes)

    @staticmethod
    def read_segmented_dataset(filename, disc_calc, segment_number):
        """Reads a segmented dataset based on specified time discontinuities."""
        try:
            with open(filename, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            print(f"Error: {filename} does not exist!")
            sys.exit()

        times = []
        magnitudes = []

        for line in lines:
            time = float(line.split()[0])
            magnitude = float(line.split()[1])

            if segment_number == 1:
                if time <= disc_calc[0]:
                    times.append(time)
                    magnitudes.append(magnitude)
            elif segment_number == len(disc_calc) + 1:
                if time >= disc_calc[-1]:
                    times.append(time)
                    magnitudes.append(magnitude)
            else:
                if disc_calc[segment_number - 2] < time <= disc_calc[segment_number - 1]:
                    times.append(time)
                    magnitudes.append(magnitude)

        print(f"There are {len(times)} events in segment {segment_number}")
        return times, magnitudes

    @staticmethod
    def adjusted_r_squared(x, y, degree):
        """Calculates the adjusted R-squared for polynomial fitting."""
        coeffs = np.polyfit(x, y, degree)
        polynomial = np.poly1d(coeffs)
        y_hat = polynomial(x)
        y_mean = np.mean(y)
        ss_reg = np.sum((y_hat - y_mean) ** 2)
        ss_tot = np.sum((y - y_mean) ** 2)
        r_squared = 1 - (((1 - (ss_reg / ss_tot)) * (len(y) - 1)) / (len(y) - degree - 1))
        return {'r_squared': r_squared}

    @staticmethod
    def read_matrix_file(filename):
        """Reads a file into a 2D numpy array (matrix)."""
        try:
            with open(filename, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            print(f"Error: {filename} does not exist!")
            sys.exit()

        matrix = np.array([[float(value) for value in line.split()] for line in lines])
        return matrix

    @staticmethod
    def read_acceptance_file(filename):
        """Reads acceptances-related metrics from a file."""
        try:
            with open(filename, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            print(f"Error: {filename} does not exist!")
            sys.exit()

        n_acc = [int(value) for value in lines[0].split()]
        n_tot = [int(value) for value in lines[1].split()]
        pct = [float(value) for value in lines[2].split()]
        return n_acc, n_tot, pct

    @staticmethod
    def read_hypercube_file(filename):
        """Reads hypercube output file with Mu, B, S, and I parameters."""
        try:
            with open(filename, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            print(f"Error: {filename} does not exist!")
            sys.exit()

        mu = [float(value) for value in lines[0].split()]
        b = [float(value) for value in lines[1].split()]
        s = [float(value) for value in lines[2].split()]
        i = [int(round(float(value))) for value in lines[3].split()]
        return mu, b, s, i