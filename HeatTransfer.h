#pragma once

#include "SimpleMath.h"

enum Example
{
	one,
};

struct HeatTransfer {

	HeatTransfer(const Example& example) {
		switch (example)
		{
		case one:
			ExampleOne();
			break;
		}
	}

	void ExampleOne() {
		double ro = 1200.0; // kg / m3
		double cp = 1500.0; // J / kgK
		double k = 0.3; // W / mK
		double L = 0.06; // m
		float Tzac = 80.0 + 273.0; // K
		float Tzrak = 20.0 + 273.0; // K
		int h = 100; // W / m2K

		int dt = 1; // s
		const int N = 6; // število vozlišè
		int t = 15000; // konèni èas v s
		double dx = L / 10;
		double diff = k / (ro * cp);
		double Fo = diff * dt / (dx * dx);
		double Bi = h * dx / k;

		sm::Vector a(N - 1);
		sm::Vector b(N);
		sm::Vector c(N - 1);
		sm::Vector z(N);

		sm::Vector T(N, Tzac);

		// define a
		for (int i = 0; i < a.m_Size; i++) {
			a[i] = -Fo;
		}
		a[a.m_Size - 1] = -2 * Fo;

		// define c
		for (int i = 0; i < c.m_Size; i++) {
			c[i] = -Fo;
		}
		c[0] = -2 * Fo;

		// define b
		for (int i = 0; i < b.m_Size; i++) {
			b[i] = 1 + 2 * Fo;
		}
		b[b.m_Size - 1] = 1 + 2 * Fo + 2 * Fo * Bi;

		int p = t / dt;
		sm::Vector temp_zg(p);
		sm::Vector temp_sp(p);
		for (int i = 0; i < p; i++) {
			for (int j = 0; j < N; j++) {
				z[j] = T[j];
			}
			z[z.m_Size - 1] = T[T.m_Size - 1] + 2 * Fo * Bi * Tzrak;

			sm::Vector Tnew = sm::algo::Thomas(a, b, c, z);
			T = Tnew;
			temp_zg[i] = T[T.m_Size - 1] - 273;
			temp_sp[i] = T[0] - 273;
		}
		PlotData(temp_sp, temp_zg);
	}
};