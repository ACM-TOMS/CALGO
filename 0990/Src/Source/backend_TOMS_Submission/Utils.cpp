/*
 This file is part of EASAL. 

 EASAL is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 EASAL is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
 * Utils.cpp
 *
 *  Created on: Jul 20, 2010
 *      Author: jrpence
 */
#include "Utils.h"
#include <fstream>
#include <Eigen/Geometry>

using namespace std;

using Eigen::Vector3f;
using Eigen::Matrix3f;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::MatrixXf;
using Eigen::VectorXf;

int Utils::findPair(vector<pair<int, int> > myvector, pair<int, int> mypair) {
	int NO = 12;
	int v1 = mypair.first;   // % (NO/2);
	int v2 = mypair.second;   // % (NO/2);

	for (int i = 0; i < myvector.size(); i++)
		if (myvector[i].first == v1 && myvector[i].second == v2)
			return i;

	return -1;
}

double Utils::dist(double *a, double *b) {
	double sum = 0.0;
	int i;
	for (i = 0; i < 3; i++) {
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}

	return sqrt(sum);
}

double Utils::dist(Vector3d a, Vector3d b) {
	double sum = 0.0;
	int i;
	for (i = 0; i < 3; i++) {
		sum += (a(i) - b(i)) * (a(i) - b(i));
	}

	return sqrt(sum);
}

double Utils::dist6(double *a, double *b) {
	double sum = 0.0;
	int i;
	for (i = 0; i < 6; i++) {
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}

	return sqrt(sum);
}

double Utils::mag(double *a) {
	double temp[3];
	temp[0] = 0;
	temp[1] = 0;
	temp[2] = 0;
	return dist(a, temp);
}

void Utils::huetoColor(double hue, float color[3]) {
	double degHue = hue * 6.0;
	int h = int(degHue) % 6;
	double f = (degHue) - floor(degHue);
	double v = 1.0;
	double p = 0;
	double q = (1 - f);
	double t = (1 - (1 - f));
	switch (h) {
	case 0:
		color[0] = v;
		color[1] = t;
		color[2] = p;
		break;
	case 1:
		color[0] = q;
		color[1] = v;
		color[2] = p;
		break;
	case 2:
		color[0] = p;
		color[1] = v;
		color[2] = t;
		break;
	case 3:
		color[0] = p;
		color[1] = q;
		color[2] = v;
		break;
	case 4:
		color[0] = t;
		color[1] = p;
		color[2] = v;
		break;
	case 5:
		color[0] = v;
		color[1] = p;
		color[2] = q;
		break;
	}

}

double* Utils::crossProd(double *b, double *c) {
	double *a = new double[3];
	a[0] = b[1] * c[2] - b[2] * c[1];
	a[1] = b[2] * c[0] - b[0] * c[2];
	a[2] = b[0] * c[1] - b[1] * c[0];
	return a;
}

double Utils::dotProd(double *a, double *b) {
	double sum = 0;
	int i;
	for (i = 0; i < 3; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

double* Utils::scalProd(double b, double* c) {
	double* a = new double[3];
	a[0] = b * c[0];
	a[1] = b * c[1];
	a[2] = b * c[2];
	return a;
}

int Utils::sign(double no) {
	if (no >= 0)
		return 1;
	else
		return -1;
}

double Utils::abs(double no) {
	if (no >= 0)
		return no;
	else
		return -no;
}

double Utils::acoss(double no) {
	if (no < -1)
		no = -1;
	else if (no > 1)
		no = 1;
	return acos(no);
}

double* Utils::vectSub(double *b, double *c) {
	double *a = new double[3];
	a[0] = b[0] - c[0];
	a[1] = b[1] - c[1];
	a[2] = b[2] - c[2];
	return a;
}

double Utils::volumeTetra(double *len) {
	double d12 = len[0] * len[0];
	double d34 = len[1] * len[1];
	double d13 = len[2] * len[2];
	double d24 = len[3] * len[3];
	double d14 = len[4] * len[4];
	double d23 = len[5] * len[5];

	//check volume constraint.
	double vc = -2 * d23 * d14 * d14
			+ (-2 * d23 * d23 + (2 * d12 + 2 * d13 + 2 * d34 + 2 * d24) * d23
					+ 2 * d24 * d13 + 2 * d12 * d34 - 2 * d12 * d24
					- 2 * d13 * d34) * d14
			+ (-2 * d34 * d24 + 2 * d24 * d13 + 2 * d12 * d34 - 2 * d12 * d13)
					* d23 - 2 * d24 * d13 * d13 + 2 * d13 * d12 * d34
			- 2 * d12 * d34 * d34 - 2 * d13 * d24 * d24 + 2 * d13 * d34 * d24
			- 2 * d12 * d12 * d34 + 2 * d12 * d34 * d24 + 2 * d12 * d24 * d13;
	return vc;
}

double Utils::volumeTetra2(double *len) {
	double u = len[0]; //d12 :0
	double U = len[1]; //d34 :1
	double v = len[2]; //d13 :2
	double V = len[3]; //d24 :3
	double w = len[4]; //d14 :4
	double W = len[5]; //d23 :5
	//U, V, W, u, v, w are lengths of edges of the tetrahedron (first three form a triangle; u opposite to U and so on),

	double X = (w - U + v) * (U + v + w);
	double x = (U - v + w) * (v - w + U);
	double Y = (u - V + w) * (V + w + u);
	double y = (V - w + u) * (w - u + V);
	double Z = (v - W + u) * (W + u + v);
	double z = (W - u + v) * (u - v + W);

	double a = sqrt((x * Y * Z));
	double b = sqrt((y * Z * X));
	double c = sqrt((z * X * Y));
	double d = sqrt((x * y * z));

	double volume = sqrt(
			((-a + b + c + d) * (a - b + c + d) * (a + b - c + d)
					* (a + b + c - d))) / (192 * u * v * w);
	return volume;
}

double dis(double *a, double *b) {
	double res = 0;
	for (int i = 0; i < 3; i++) {
		res += (b[i] - a[i]) * (b[i] - a[i]);
	}
	return sqrt(res);
}

double* Utils::newMatApp(double* x, double* a, double* b, double* c, double* ap,
		double* bp, double* cp) {
	double x1 = x[0];
	double x2 = x[1];
	double x3 = x[2];

	double a1 = a[0];
	double a2 = a[1];
	double a3 = a[2];

	double b1 = b[0];
	double b2 = b[1];
	double b3 = b[2];

	double c1 = c[0];
	double c2 = c[1];
	double c3 = c[2];

	double ap1 = ap[0];
	double ap2 = ap[1];
	double ap3 = ap[2];

	double bp1 = bp[0];
	double bp2 = bp[1];
	double bp3 = bp[2];

	double cp1 = cp[0];
	double cp2 = cp[1];
	double cp3 = cp[2];

	//the output
	double* cg3 = new double[3];
	double t1 = x1 - a1;
	double t2 = bp1 - ap1;
	double t4 = x2 - a2;
	double t5 = bp2 - ap2;
	double t7 = x3 - a3;
	double t8 = bp3 - ap3;
	double t10 = t1 * t2 + t4 * t5 + t7 * t8;
	double t11 = c2 * c2;
	double t13 = b1 * c2;
	double t14 = t13 * a2;
	double t17 = c2 * b2;
	double t19 = t17 * a1;
	double t20 = a2 * c2;
	double t22 = a2 * a2;
	double t23 = b1 * t22;
	double t24 = a2 * a1;
	double t25 = t24 * c2;
	double t26 = a2 * b2;
	double t27 = t26 * c1;
	double t28 = t26 * a1;
	double t29 = t22 * c1;
	double t30 = -b1 * t11 + 0.2e1 * t14 + a1 * t11 + t17 * c1 - t19 - t20 * c1
			- t23 - t25 - t27 + t28 + t29;
	double t31 = b3 * c1;
	double t33 = t31 * a3;
	double t34 = b3 * a1;
	double t35 = t34 * c3;
	double t36 = t34 * a3;
	double t37 = a3 * c1;
	double t39 = a3 * a3;
	double t40 = t39 * c1;
	double t41 = c3 * c3;
	double t43 = b1 * c3;
	double t44 = t43 * a3;
	double t46 = b1 * t39;
	double t48 = a1 * c3;
	double t49 = t48 * a3;
	double t50 = t31 * c3 - t33 - t35 + t36 - t37 * c3 + t40 - b1 * t41
			+ 0.2e1 * t44 - t46 + a1 * t41 - t49;
	double t54 = b1 * a3;
	double t60 = c1 * c1;
	double t61 = b2 * t60;
	double t63 = b3 * b3;
	double t64 = t63 * c1;
	double t67 = b2 * b2;
	double t68 = t67 * c3;
	double t71 = a1 * b1;
	double t73 = a1 * a1;
	double t74 = b2 * t73;
	double t76 = t67 * c1;
	double t78 = b3 * t60;
	double t80 = -t31 * t43 + t31 * t54 + t31 * t48 + t34 * t37 + t34 * t43
			- t34 * t54 - t61 * a2 - t64 * a1 - t29 * b1 - t68 * a3 - t26 * t41
			- t71 * t11 - t74 * c2 - t76 * a1 - t78 * a3;
	double t81 = b3 * t73;
	double t86 = b3 * t22;
	double t89 = t63 * c2;
	double t92 = b3 * a3;
	double t95 = b2 * t39;
	double t98 = b1 * b1;
	double t99 = c2 * t98;
	double t102 = t98 * c3;
	double t114 = -0.2e1 * t81 * c3 - 0.2e1 * t40 * b1 - 0.2e1 * t86 * c3
			- 0.2e1 * t89 * a2 - 0.2e1 * t92 * t11 - 0.2e1 * t95 * c2
			- 0.2e1 * t99 * a2 - 0.2e1 * t102 * a3 - 0.2e1 * t71 * t41
			+ t63 * t60 + t63 * t73 + t39 * t60 + t63 * t11 + t63 * t22
			+ t39 * t11 + t67 * t60;
	double t131 = a3 * a1;
	double t134 = b3 * b2;
	double t135 = c3 * c2;
	double t138 = a2 * c3;
	double t141 = t67 * t73 + t22 * t60 + t67 * t41 + t67 * t39 + t22 * t41
			+ t98 * t11 + t98 * t22 + t98 * t41 + t98 * t39 + t73 * t11
			+ t73 * t41 + 0.2e1 * t37 * t43 - 0.2e1 * t37 * t48
			+ 0.2e1 * t131 * t43 - 0.2e1 * t134 * t135 + 0.2e1 * t134 * t138;
	double t142 = a3 * c2;
	double t144 = a2 * a3;
	double t146 = b3 * a2;
	double t149 = b2 * a3;
	double t153 = b2 * c1;
	double t155 = b1 * a2;
	double t157 = a1 * c2;
	double t159 = b2 * a1;
	double t162 = a2 * c1;
	double t167 = t134 * t142 - t134 * t144 + t146 * t135 + t92 * t20
			+ t149 * t135 + t149 * t138 - t144 * t135 - t153 * t13 + t153 * t155
			+ t153 * t157 + t159 * t13 - t159 * t155 + t159 * t162 + t162 * t13
			- t162 * t157 + t24 * t13;
	double t170 = 0.1e1 / (0.2e1 * t80 + t114 + t141 + 0.2e1 * t167);
	double t172 = cp1 - ap1;
	double t174 = cp2 - ap2;
	double t176 = cp3 - ap3;
	double t178 = t1 * t172 + t4 * t174 + t7 * t176;
	double t179 = b2 * b1;
	double t184 = -t179 * c2 + t179 * a2 + t19 + t76 - t67 * a1 - 0.2e1 * t27
			+ t14 - t23 - t25 + t28 + t29;
	double t187 = b3 * b1;
	double t190 = t64 - 0.2e1 * t33 - t63 * a1 + t36 + t40 - t187 * c3 + t44
			+ t187 * a3 - t46 + t35 - t49;
	double t206 = t1 * (t5 * t176 - t8 * t174) + t4 * (t8 * t172 - t2 * t176)
			+ t7 * (t2 * t174 - t5 * t172);
	double t207 = b2 * c3;
	double t213 = c1 * b1;
	double t215 = t213 * a2;
	double t216 = a1 * c1;
	double t218 = t153 * a1;
	double t221 = t71 * c2;
	double t222 = t71 * a2;
	double t223 = t73 * c2;
	double t224 = t24 * c1;
	double t225 = -t213 * c2 + t215 + t216 * c2 + t61 - 0.2e1 * t218 - a2 * t60
			+ t221 - t222 - t223 + t74 + t224;
	double t227 = t149 * c3;
	double t230 = t144 * c3;
	double t233 = t92 * c2;
	double t234 = t146 * c3;
	double t235 = t92 * a2;
	double t238 = t39 * c2;
	double t239 = b2 * t41 - 0.2e1 * t227 + t95 - a2 * t41 + t230 - b3 * c3 * c2
			+ t233 + t234 - t235 + a3 * c3 * c2 - t238;
	double t247 = -t99 + t98 * a2 + 0.2e1 * t221 + t153 * b1 - t179 * a1 - t215
			- t222 - t223 - t218 + t74 + t224;
	double t252 = t134 * c3 - t134 * a3 - t234 - t89 + t63 * a2 + 0.2e1 * t233
			- t227 + t95 + t230 - t235 - t238;
	double t260 = t31 * a1;
	double t264 = t213 * a3;
	double t266 = t131 * c1;
	double t267 = t71 * c3;
	double t268 = t71 * a3;
	double t269 = t73 * c3;
	double t270 = t78 - 0.2e1 * t260 - a3 * t60 - t213 * c3 + t264 + t216 * c3
			+ t81 + t266 + t267 - t268 - t269;
	double t272 = t26 * c3;
	double t273 = t149 * c2;
	double t274 = t26 * a3;
	double t276 = t22 * c3;
	double t278 = t146 * c2;
	double t281 = t144 * c2;
	double t282 = -t207 * c2 + t272 + t273 - t274 + t138 * c2 - t276 + b3 * t11
			- 0.2e1 * t278 + t86 - a3 * t11 + t281;
	double t290 = t187 * c1 - t187 * a1 - t264 - t102 + t98 * a3 + 0.2e1 * t267
			- t260 + t81 + t266 - t268 - t269;
	double t295 = -t68 + t67 * a3 + 0.2e1 * t272 + t134 * c2 - t134 * a2 - t273
			- t274 - t276 - t278 + t86 + t281;
	cg3[0] = ap1 + -t10 * (t30 + t50) * t170 + t178 * (t184 + t190) * t170
			+ t206 * (t207 - t149 - t138 - b3 * c2 + t146 + t142) * t170;
	cg3[1] = ap2 + t10 * (t225 + t239) * t170 - t178 * (t247 + t252) * t170
			- t206 * (-t31 + t34 + t37 + t43 - t54 - t48) * t170;
	cg3[2] = ap3 + t10 * (t270 + t282) * t170 - t178 * (t290 + t295) * t170
			+ t206 * (t13 - t155 - t157 - t153 + t159 + t162) * t170;

	return cg3;
}

double* Utils::find4thPoint_inTheSamePlane(double* pa, double* pb, double* pc,
		double ax, double bx, double cx) {
	Vector3d result(0, 0, 0);
	Vector3d zero(0, 0, 0);
	Vector3d a(pa[0], pa[1], pa[2]);
	Vector3d b(pb[0], pb[1], pb[2]);
	Vector3d c(pc[0], pc[1], pc[2]);

	double ab = Utils::dist(pa, pb);
	double bc = Utils::dist(pb, pc);
	double ac = Utils::dist(pa, pc);

	double cosXAB = (ax * ax + ab * ab - bx * bx) / (2 * ax * ab);
	double cosXAC = (ax * ax + ac * ac - cx * cx) / (2 * ax * ac);
	double cosXBC = (bx * bx + bc * bc - cx * cx) / (2 * bx * bc);

	//if the chosen triangle is kind of line then calculations turn wrong, so chose the triangle as possible not to be a line.
	double min = abs(abs(cosXAB) - 0.5); //the cos that is closest to middle angle
	int mcase = 1;
	if (min > abs(abs(cosXAC) - 0.5)) {
		min = abs(abs(cosXAC) - 0.5);
		mcase = 2;
	}
	if (min > abs(abs(cosXBC) - 0.5)) {
		min = abs(abs(cosXBC) - 0.5);
		mcase = 3;
	}
	if (mcase == 1)  //cosXAB
			{
		double sinXAB = sqrt(abs(1 - cosXAB * cosXAB));

		Vector3d vab = b - a;
		Vector3d vac = c - a;
		Vector3d v3 = vab.cross(vac); //perpendicular to abc plane
		//c = a X b : the first (index) finger can represent a, the first vector in the product; the second (middle) finger, b, the second vector; and the thumb, c, the product
		Vector3d vab_perp = v3.cross(vab); //perpendicular to vab vector
		double vab_p = Utils::dist(vab_perp, zero);

		Vector3d r = a + ax * cosXAB * vab / ab; //walk along vab direction first ax*cosXAB amount
		Vector3d r1 = r + ax * sinXAB * vab_perp / vab_p; //walk along perpendicular to vab direction ax*sinXAB amount
		Vector3d r2 = r - ax * sinXAB * vab_perp / vab_p; //walk along opposite perpendicular to vab direction ax*sinXAB amount

		Vector3d r1c = r1 - c;
		double dr1c = Utils::dist(r1c, zero);

		Vector3d r2c = r2 - c;
		double dr2c = Utils::dist(r2c, zero);
		if (abs(dr1c - cx) > abs(dr2c - cx)) //the other side
			result = r2;
		else
			result = r1;
	} else if (mcase == 2)  //cosXAC
			{
		double sinXAC = sqrt(abs(1 - cosXAC * cosXAC));

		Vector3d vab = b - a;
		Vector3d vac = c - a;
		Vector3d v3 = vac.cross(vab); //perpendicular to abc plane
		//c = a X b : the first (index) finger can represent a, the first vector in the product; the second (middle) finger, b, the second vector; and the thumb, c, the product
		Vector3d vac_perp = v3.cross(vac); //perpendicular to vab vector
		double vac_p = Utils::dist(vac_perp, zero);

		Vector3d r = a + ax * cosXAC * vac / ac; //walk along vac direction first ax*cosXAC amount
		Vector3d r1 = r + ax * sinXAC * vac_perp / vac_p; //walk along perpendicular to vac direction ax*sinXAC amount
		Vector3d r2 = r - ax * sinXAC * vac_perp / vac_p; //walk along opposite perpendicular to vac direction ax*sinXAC amount

		Vector3d r1b = r1 - b;
		double dr1b = Utils::dist(r1b, zero);

		Vector3d r2b = r2 - b;
		double dr2b = Utils::dist(r2b, zero);
		if (abs(dr1b - bx) > abs(dr2b - bx)) //the other side
			result = r2;
		else
			result = r1;
	} else if (mcase == 3)  //cosXBC
			{
		double sinXBC = sqrt(abs(1 - cosXBC * cosXBC));

		Vector3d vbc = c - b;
		Vector3d vba = a - b;
		Vector3d v3 = vbc.cross(vba); //perpendicular to abc plane
		//c = a X b : the first (index) finger can represent a, the first vector in the product; the second (middle) finger, b, the second vector; and the thumb, c, the product
		Vector3d vbc_perp = v3.cross(vbc); //perpendicular to vbc vector
		double vbc_p = Utils::dist(vbc_perp, zero);

		Vector3d r = b + bx * cosXBC * vbc / bc; //walk along vbc direction first bx*cosXBC amount
		Vector3d r1 = r + bx * sinXBC * vbc_perp / vbc_p; //walk along perpendicular to vbc direction bx*sinXBC amount
		Vector3d r2 = r - bx * sinXBC * vbc_perp / vbc_p; //walk along opposite perpendicular to vbc direction bx*sinXBC amount

		Vector3d r1a = r1 - a;
		double dr1a = Utils::dist(r1a, zero);

		Vector3d r2a = r2 - a;
		double dr2a = Utils::dist(r2a, zero);
		if (abs(dr1a - ax) > abs(dr2a - ax)) //the other side
			result = r2;
		else
			result = r1;
	}

	double* cg3 = new double[3];
	cg3[0] = result(0);
	cg3[1] = result(1);
	cg3[2] = result(2);
	return cg3;
}

double* Utils::matApp2(double* x, double* a, double* b, double* c, double* ap,
		double* bp, double* cp) {

	Vector3d p1(a[0], a[1], a[2]);
	Vector3d p2(b[0], b[1], b[2]);
	Vector3d p3(c[0], c[1], c[2]);

	Vector3d P1(ap[0], ap[1], ap[2]);
	Vector3d P2(bp[0], bp[1], bp[2]);
	Vector3d P3(cp[0], cp[1], cp[2]);

	Vector3d v1 = p2 - p1;
	Vector3d v2 = p3 - p1;
	Vector3d v3 = v1.cross(v2);

	Vector3d V1 = P2 - P1;
	Vector3d V2 = P3 - P1;
	Vector3d V3 = V1.cross(V2);

	Matrix3d m, M, R;
	m << v1(0), v2(0), v3(0), v1(1), v2(1), v3(1), v1(2), v2(2), v3(2);
	M << V1(0), V2(0), V3(0), V1(1), V2(1), V3(1), V1(2), V2(2), V3(2);

	R = M * m.inverse(); //rotation matrix
	Vector3d t = P1 - R * p1;  //translation

	Vector3d p(x[0], x[1], x[2]);
	Vector3d CG3 = R * p + t;

	double* cg3 = new double[3];
	cg3[0] = CG3(0);
	cg3[1] = CG3(1);
	cg3[2] = CG3(2);
	return cg3;
}

Vector3d Utils::RotMatToEuler(Matrix3d rmat) {

	double phi = 0.0;
	double cos_theta = 1.0;
	double psi = 0.0;

	cos_theta = rmat(2, 2);

	if (abs(cos_theta - 1.0) > 1.0e-8) {
		psi = atan2(rmat(0, 2), rmat(1, 2));
		phi = atan2(rmat(2, 0), -rmat(2, 1));
	} else {
		psi = atan2(-rmat(1, 0), rmat(0, 0));
		phi = 0.0;
	}

	Vector3d out(phi, cos_theta, psi);
	return out;
}

void Utils::matApp(double* x, double* a, double* b, double* c, double* ap,
		double* bp, double* cp, double cg3[3]) {
	double x1 = x[0];
	double x2 = x[1];
	double x3 = x[2];

	double a1 = a[0];
	double a2 = a[1];
	double a3 = a[2];

	double b1 = b[0];
	double b2 = b[1];
	double b3 = b[2];

	double c1 = c[0];
	double c2 = c[1];
	double c3 = c[2];

	double ap1 = ap[0];
	double ap2 = ap[1];
	double ap3 = ap[2];

	double bp1 = bp[0];
	double bp2 = bp[1];
	double bp3 = bp[2];

	double cp1 = cp[0];
	double cp2 = cp[1];
	double cp3 = cp[2];

	//the output
	double t1 = (b1 * b1);
	double t2 = (ap1 * t1);
	double t3 = (c3 * c3);
	double t5 = (b2 * b2);
	double t6 = (ap1 * t5);
	double t7 = (c1 * c1);
	double t9 = (b3 * b3);
	double t10 = (ap1 * t9);
	double t11 = (c2 * c2);
	double t13 = (a3 * a3);
	double t14 = (bp1 * t13);
	double t17 = (cp1 * t9);
	double t18 = (a2 * a2);
	double t20 = (a1 * a1);
	double t21 = (bp1 * t20);
	double t23 = (bp1 * b3);
	double t24 = (a2 * c3);
	double t25 = (t24 * x2);
	double t27 = (ap1 * b3);
	double t28 = (c1 * a1);
	double t29 = (t28 * c3);
	double t31 = (ap3 * cp2);
	double t32 = (b1 * c2);
	double t33 = (t32 * a3);
	double t35 = (ap1 * c2);
	double t36 = (a2 * c1);
	double t37 = (t36 * x1);
	double t39 = (bp1 * a1);
	double t40 = (b1 * a3);
	double t41 = (t40 * x3);
	double t43 = (ap2 * cp3);
	double t44 = (b3 * a2);
	double t45 = (t44 * x1);
	double t47 = (-t2 * t3 - t6 * t7 - t10 * t11 - t14 * t11 - t2 * t11
			- t17 * t18 - t21 * t3 - t23 * t25 - t27 * t29 + t31 * t33
			+ t35 * t37 + t39 * t41 + t43 * t45);
	double t48 = (b2 * a1);
	double t49 = (t48 * x3);
	double t51 = (t36 * x3);
	double t53 = (t32 * x2);
	double t55 = (bp3 * cp2);
	double t56 = (b2 * c1);
	double t57 = (t56 * a3);
	double t59 = (b3 * c1);
	double t60 = (t59 * x2);
	double t62 = (a3 * c2);
	double t63 = (t62 * x1);
	double t65 = b3 * c2;
	double t66 = (t65 * a1);
	double t68 = (bp2 * ap3);
	double t71 = (c3 * c2);
	double t72 = (t71 * x2);
	double t74 = (b2 * c3);
	double t75 = (t74 * x1);
	double t78 = (t32 * x3);
	double t80 = bp3 * ap2;
	double t82 = (t43 * t49 + t43 * t51 - t39 * t53 + t55 * t57 + t55 * t60
			+ t43 * t63 + t43 * t66 + t68 * t66 + t68 * t45 + t23 * t72
			+ t55 * t75 + t55 * t51 + t55 * t78 - t80 * t49);
	double t85 = (t56 * x3);
	double t87 = (a1 * c2);
	double t88 = (t87 * x3);
	double t90 = (b1 * a2);
	double t91 = (t90 * x3);
	double t93 = (t40 * x2);
	double t96 = (a1 * c3);
	double t97 = (t96 * x2);
	double t100 = (b1 * c3);
	double t101 = (t100 * a2);
	double t103 = (t100 * x2);
	double t105 = (a3 * c1);
	double t106 = (t105 * x2);
	double t108 = (b3 * a1);
	double t109 = (t108 * x2);
	double t111 = (t59 * a2);
	double t114 = (-t80 * t57 + t80 * t85 + t80 * t88 + t80 * t91 + t68 * t93
			+ t80 * t33 - t80 * t97 - t80 * t93 - t80 * t101 + t80 * t103
			+ t80 * t106 + t80 * t109 + t80 * t111 - t80 * t60);
	double t118 = (t65 * x1);
	double t120 = (t24 * x1);
	double t122 = (b2 * a3);
	double t123 = (t122 * x1);
	double t125 = (t74 * a1);
	double t134 = (-t80 * t63 - t80 * t45 - t80 * t66 + t80 * t118 + t80 * t120
			+ t80 * t123 + t80 * t125 - t80 * t75 - t80 * t51 - t80 * t78
			- t31 * t49 - t31 * t57 + t31 * t85 + t31 * t88);
	double t150 = (t31 * t91 - t31 * t97 - t31 * t93 - t31 * t101 + t31 * t103
			+ t31 * t106 + t31 * t109 + t31 * t111 - t31 * t60 - t31 * t63
			- t31 * t45 - t31 * t66 + t31 * t118);
	double t157 = (t74 * c2);
	double t160 = (ap1 * b2);
	double t161 = (c1 * b1);
	double t162 = (t161 * c2);
	double t165 = (cp1 * b3);
	double t166 = (t161 * a3);
	double t168 = (cp1 * a3);
	double t169 = (a1 * b1);
	double t170 = (t169 * c3);
	double t172 = (t169 * a3);
	double t176 = (a1 * a3);
	double t177 = (t176 * c1);
	double t179 = (t74 * a2);
	double t181 = t31 * t120 + t31 * t123 + t31 * t125 - t31 * t75 - t31 * t51
			- t31 * t78 + 2 * t27 * t157 + 2 * t160 * t162 - t165 * t166
			- t168 * t170 + 2 * t165 * t172 - t165 * t170 - t165 * t177
			- t165 * t179;
	double t183 = (t62 * a2);
	double t185 = (t122 * a2);
	double t188 = (t122 * c2);
	double t190 = (cp1 * c2);
	double t191 = (t90 * x1);
	double t193 = (ap1 * a2);
	double t194 = (t56 * x1);
	double t196 = (t28 * c2);
	double t199 = (c1 * c3 * x1);
	double t201 = (t96 * x1);
	double t205 = (t105 * x1);
	double t207 = (ap1 * a3);
	double t209 = (ap1 * b1);
	double t210 = (c3 * a3);
	double t211 = (t210 * x1);
	double t215 = -t165 * t183 + 2 * t165 * t185 - t165 * t188 - t190 * t191
			- t193 * t194 - t160 * t196 - t27 * t199 + 2 * t27 * t201
			- t193 * t162 - t27 * t205 + t207 * t199 - t209 * t211 + t43 * t75
			+ t43 * t78;
	double t230 = (t55 * t49 - t55 * t85 - t55 * t88 - t55 * t91 - t55 * t33
			+ t55 * t97 + t55 * t93 + t55 * t101 - t55 * t103 - t55 * t106
			- t55 * t109 - t55 * t111 + t55 * t63 + t55 * t45);
	double t240 = (cp1 * b2);
	double t241 = (t161 * a2);
	double t243 = (t161 * c3);
	double t247 = (bp1 * t18);
	double t249 = (cp1 * t1);
	double t253 = t55 * t66 - t55 * t118 - t55 * t120 - t55 * t123 - t55 * t125
			- t23 * t201 - t240 * t241 + 2 * t27 * t243 + t207 * t72 - t247 * t3
			- t249 * t13 - t21 * t11 - t14 * t7;
	double t256 = (cp1 * t5);
	double t259 = (bp1 * b2);
	double t260 = t3 * x2;
	double t262 = (bp1 * c1);
	double t263 = (t96 * x3);
	double t265 = (t105 * x3);
	double t267 = (t100 * x3);
	double t269 = t48 * x2;
	double t272 = t87 * x2;
	double t274 = (t71 * x3);
	double t276 = (t24 * x3);
	double t278 = (t62 * x3);
	double t280 = (bp1 * a2);
	double t281 = (t122 * x3);
	double t283 = -t6 * t3 - t247 * t7 - t256 * t13 - t256 * t20 - t259 * t260
			- t262 * t263 - t39 * t265 - t39 * t267 + 2 * t262 * t269
			- t262 * t272 + t259 * t274 - t259 * t276 - t259 * t278
			+ t280 * t281;
	double t286 = (t65 * x3);
	double t290 = (bp1 * a3);
	double t292 = t74 * x2;
	double t295 = t36 * x2;
	double t297 = t90 * x2;
	double t299 = (t108 * x3);
	double t303 = t44 * x2;
	double t305 = t62 * x2;
	double t309 = (t24 * c2);
	double t311 = -t280 * t274 + 2 * t280 * t286 - t280 * t278 - t290 * t25
			+ 2 * t290 * t292 - t39 * t295 + t39 * t297 + 2 * t262 * t299
			- t290 * t72 + t290 * t303 - t23 * t305 + t262 * t267 - t207 * t243
			- t27 * t309;
	double t314 = t48 * x1;
	double t328 = -t35 * t191 - t35 * t194 + 2 * t35 * t314 - t207 * t157
			- t43 * t103 - t43 * t106 - t43 * t109 - t43 * t111 + t43 * t60
			- t43 * t118 - t43 * t120 - t43 * t123 - t43 * t125 - t23 * t29;
	double t338 = (bp1 * c2);
	double t346 = t87 * x1;
	double t348 = -t290 * t170 + 2 * t290 * t29 - t290 * t243 - t23 * t177
			- t23 * t183 - t23 * t309 + 2 * t338 * t191 + t338 * t194
			- t338 * t37 - t338 * t314 - t290 * t157 - t280 * t194
			- t280 * t346;
	double t356 = (t169 * c2);
	double t359 = (ap1 * c1);
	double t366 = (ap1 * a1);
	double t370 = 2 * t290 * t309 - t290 * t179 - t259 * t196 + 2 * t280 * t196
			+ t23 * t199 - t280 * t356 + t280 * t314 + 2 * t359 * t297
			- t359 * t53 + 2 * t359 * t41 + t359 * t263 - t366 * t267
			- t359 * t269 + t359 * t272;
	double t387 = (a1 * a2);
	double t388 = (t387 * c1);
	double t390 = -t366 * t53 - t160 * t274 - t160 * t276 + 2 * t160 * t278
			+ t193 * t274 - t193 * t286 - t27 * t72 - t207 * t292 - t359 * t299
			+ 2 * t27 * t25 - t27 * t305 - t359 * t267 - t280 * t162
			- t259 * t388;
	double t393 = t176 * x1;
	double t396 = (bp1 * b1);
	double t403 = (cp1 * a2);
	double t410 = -t23 * t205 - t290 * t199 + t23 * t393 - t39 * t211
			+ 2 * t396 * t211 - t262 * t297 + t262 * t53 - t262 * t41
			- t190 * t314 + 2 * t403 * t194 + t403 * t346 - t168 * t179
			- t240 * t356 + t240 * t286;
	double t415 = t44 * x3;
	double t419 = (t169 * a2);
	double t426 = t122 * x2;
	double t428 = (cp1 * b1);
	double t429 = t59 * x3;
	double t433 = (cp1 * a1);
	double t436 = -t240 * t415 - t403 * t356 - t403 * t314 + 2 * t240 * t419
			- t165 * t201 - t240 * t388 + 2 * t165 * t205 - t165 * t426
			+ t428 * t429 - t428 * t299 - t165 * t393 + t433 * t211
			- t428 * t211;
	double t439 = t32 * x1;
	double t442 = t100 * x1;
	double t444 = (cp1 * c1);
	double t451 = t40 * x1;
	double t453 = t161 * x2;
	double t458 = -t428 * t269 + t165 * t292 + t240 * t439 - t240 * t191
			+ t165 * t442 - t444 * t297 - t444 * t41 + t433 * t265
			+ 2 * t433 * t267 - t444 * t269 - t165 * t451 + t240 * t453
			+ 2 * t433 * t53 - t433 * t41;
	double t476 = 2 * t240 * t276 - t240 * t278 - t403 * t281 - t403 * t286
			+ t403 * t278 + t168 * t25 - t168 * t292 + t433 * t295 - t433 * t297
			- t444 * t299 - t168 * t303 - t165 * t25 + 2 * t165 * t305
			- t27 * t166;
	double t481 = c2 * x2;
	double t483 = t20 * x3;
	double t485 = t7 * a3;
	double t487 = a3 * t11;
	double t489 = a2 * t3;
	double t491 = a1 * t11;
	double t493 = t7 * a2;
	double t495 = a1 * t3;
	double t497 = t20 * x2;
	double t499 = cp1 * t20;
	double t501 = -t27 * t170 - t27 * t179 - t27 * t188 + t6 * t210 + t10 * t481
			- t23 * t483 + t27 * t485 + t27 * t487 + t160 * t489 + t209 * t491
			+ t160 * t493 + t209 * t495 + t240 * t497 - t499 * t481;
	double t504 = t13 * x1;
	double t506 = cp1 * t13;
	double t507 = c1 * x1;
	double t509 = cp1 * t18;
	double t511 = t18 * x1;
	double t513 = c2 * a2;
	double t515 = t18 * x3;
	double t517 = c3 * x3;
	double t521 = t13 * x2;
	double t523 = a2 * x2;
	double t526 = t1 * x2;
	double t528 = t428 * t504 - t506 * t507 - t509 * t507 + t428 * t511
			+ t17 * t513 + t165 * t515 - t509 * t517 - t499 * t517 - t506 * t481
			+ t240 * t521 + t17 * t523 + t249 * t523 - t190 * t526;
	double t529 = a1 * x1;
	double t535 = t20 * c3;
	double t538 = a3 * x3;
	double t540 = t18 * c3;
	double t546 = t13 * c2;
	double t548 = t17 * t529 - t17 * t507 + t256 * t529 - t256 * t507
			+ t506 * t161 + t165 * t535 + t17 * t28 + t256 * t538 + t165 * t540
			- t256 * t517 + t249 * t538 - t249 * t517 + t256 * t28
			+ t240 * t546;
	double t551 = t20 * c2;
	double t554 = t1 * a2;
	double t566 = t509 * t161 + t240 * t551 + t249 * t210 + t190 * t554
			+ t256 * t210 - t17 * t481 + t165 * t483 + t10 * t513 - t10 * t523
			- t2 * t523 + t35 * t526 - t10 * t529 + t10 * t507 - t6 * t529;
	double t577 = t7 * x2;
	double t580 = t3 * x1;
	double t583 = t6 * t507 + t10 * t28 - t6 * t538 + t6 * t517 - t2 * t538
			+ t2 * t517 + t6 * t28 + t2 * t210 - t259 * t497 + t21 * t481
			+ t280 * t577 - t259 * t577 + t39 * t580 - t396 * t504;
	double t591 = t11 * x1;
	double t594 = t11 * x3;
	double t600 = t7 * x3;
	double t603 = -t396 * t580 + t14 * t507 + t247 * t507 - t396 * t511
			+ t39 * t591 - t396 * t591 + t290 * t594 - t23 * t515 - t23 * t594
			+ t247 * t517 + t21 * t517 + t290 * t600 - t23 * t600;
	double t618 = t14 * t481 + t280 * t260 - t259 * t521 + t14 * t161
			+ t23 * t535 + t23 * t485 + t23 * t487 + t23 * t540 + t259 * t489
			+ t259 * t546 + t396 * t491 + t247 * t161 + t259 * t551
			+ t259 * t493;
	double t634 = t396 * t495 + t160 * t260 - t193 * t577 + t160 * t577
			- t366 * t580 + t209 * t580 - t366 * t591 + t209 * t591
			- t207 * t594 + t27 * t594 - t207 * t600 + t27 * t600 - t193 * t260
			+ t35 * t554;
	double t649 = -t10 * t7 - t17 * t20 - t249 * t18 - t160 * t356 - t160 * t241
			- t160 * t286 + t160 * t415 + t27 * t426 - t209 * t429 + t209 * t299
			+ t209 * t269 - t27 * t292 - t160 * t439 + t160 * t191;
	double t655 = (bp2 * cp3);
	double t667 = -t27 * t442 + t27 * t451 - t160 * t453 - t655 * t49
			- t655 * t57 + t655 * t85 + t655 * t88 + t655 * t91 + t655 * t33
			- t655 * t97 - t655 * t93 - t655 * t101 + t655 * t103 + t655 * t106;
	double t682 = t655 * t109 + t655 * t111 - t655 * t60 - t655 * t63
			- t655 * t45 - t655 * t66 + t655 * t118 + t655 * t120 + t655 * t123
			+ t655 * t125 - t655 * t75 - t655 * t51 - t655 * t78 + t68 * t49;
	double t698 = t68 * t57 - t68 * t85 - t68 * t88 - t68 * t91 - t68 * t33
			+ t68 * t97 + t68 * t101 - t68 * t103 - t68 * t106 - t68 * t109
			- t68 * t111 + t68 * t60 + t68 * t63 - t68 * t118;
	double t713 = -t68 * t120 - t68 * t123 - t68 * t125 + t68 * t75 + t68 * t51
			+ t68 * t78 + t43 * t57 - t43 * t85 - t43 * t88 - t43 * t91
			- t43 * t33 + t43 * t97 + t43 * t93 + t43 * t101;
	double t748 = -2 * t5 * c1 * a1 - 2 * t18 * c1 * b1 - 2 * t13 * c1 * b1
			- 2 * c2 * t1 * a2 - 2 * b2 * t7 * a2 - 2 * b2 * t13 * c2 + t5 * t3
			+ t18 * t3 + t20 * t11 - 2 * t1 * c3 * a3 + t1 * t11 + t13 * t7
			+ t18 * t7 + t1 * t3 + t5 * t20;
	double t776 = t9 * t11 + t5 * t13 + t1 * t18 - 2 * b3 * t7 * a3 + t5 * t7
			+ t13 * t11 + t1 * t13 - 2 * t59 * t100 + 2 * t59 * t96
			+ 2 * t59 * t40 + 2 * t176 * t100 - 2 * t105 * t96 + 2 * t105 * t100
			- 2 * t108 * t40 + 2 * t108 * t100 + 2 * t108 * t105;
	double t778 = b3 * b2;
	double t781 = b3 * a3;
	double t784 = a3 * a2;
	double t797 = t778 * t24 - t778 * t71 + t781 * t513 + t44 * t71
			- t778 * t784 + t778 * t62 + t122 * t71 - t784 * t71 + t122 * t24
			+ t48 * t32 + t56 * t87 + t56 * t90 - t56 * t32 - t36 * t87
			+ t387 * t32 - t48 * t90;
	double t833 = 2 * t36 * t32 + 2 * t48 * t36 + t9 * t18 + t9 * t20 + t9 * t7
			- 2 * t9 * c2 * a2 - 2 * t9 * c1 * a1 - 2 * t5 * c3 * a3
			- 2 * t781 * t11 - 2 * b3 * t18 * c3 - 2 * b2 * t20 * c2
			- 2 * b3 * t20 * c3 - 2 * t169 * t11 - 2 * b2 * a2 * t3 + t20 * t3
			- 2 * t169 * t3;
	double t836 = 1 / (t748 + t776 + 2 * t797 + t833);
	double t838 = cp2 * t5;
	double t840 = cp2 * t1;
	double t842 = ap2 * t5;
	double t844 = bp2 * t20;
	double t846 = bp2 * t13;
	double t848 = cp2 * a2;
	double t851 = (ap1 * cp3);
	double t853 = bp3 * cp1;
	double t855 = bp2 * a2;
	double t858 = bp3 * ap1;
	double t860 = ap3 * cp1;
	double t862 = -t838 * t13 - t840 * t18 - t842 * t3 - t844 * t11 - t846 * t11
			- t848 * t281 - t846 * t7 + t851 * t111 + t853 * t120 - t855 * t346
			+ t848 * t346 + t858 * t97 + t860 * t60;
	double t866 = ap2 * b1;
	double t870 = bp2 * a1;
	double t872 = cp2 * b2;
	double t874 = bp2 * a3;
	double t877 = bp2 * b3;
	double t882 = t860 * t51 + t851 * t103 + t853 * t118 + t866 * t299
			+ t858 * t60 + t858 * t66 + t870 * t41 + t872 * t453 + t874 * t303
			+ t860 * t66 - t877 * t29 + t853 * t109 + t877 * t199 + t853 * t125;
	double t898 = -t853 * t75 - t853 * t51 - t853 * t78 + t858 * t49
			+ t858 * t57 - t858 * t85 - t858 * t88 - t858 * t91 - t858 * t33
			+ t858 * t93 - t858 * t103 - t858 * t106 - t858 * t109
			- t858 * t111;
	double t905 = ap2 * b3;
	double t910 = ap2 * b2;
	double t916 = t858 * t45 - t874 * t170 + 2 * t874 * t29 - t874 * t243
			- t877 * t177 + t905 * t426 - t866 * t429 + t866 * t269
			- t905 * t292 - t910 * t439 + t910 * t191 - t905 * t442
			- t905 * t166 - t905 * t170;
	double t932 = -t905 * t179 - t905 * t188 - t910 * t356 - t910 * t241
			- t910 * t286 + t910 * t415 - t858 * t118 - t858 * t120
			- t858 * t123 - t858 * t125 + t858 * t75 + t858 * t51 + t858 * t78;
	double t945 = bp1 * ap3;
	double t948 = t860 * t49 + t860 * t57 - t860 * t85 - t860 * t88 - t860 * t91
			- t860 * t33 + t860 * t97 + t860 * t93 + t860 * t101 - t860 * t103
			- t860 * t106 - t860 * t109 + t945 * t33 - t945 * t97;
	double t964 = -t945 * t93 - t945 * t101 + t945 * t103 + t945 * t106
			+ t945 * t109 + t945 * t111 - t945 * t60 - t945 * t63 + t853 * t111
			+ t853 * t106 + t851 * t120 - t877 * t309 - t851 * t60
			+ t945 * t118;
	double t965 = ap2 * t9;
	double t967 = ap2 * t1;
	double t969 = ap2 * c2;
	double t978 = cp2 * t20;
	double t980 = cp2 * b1;
	double t982 = cp2 * t13;
	double t984 = cp2 * t18;
	double t986 = -t965 * t523 - t967 * t523 + t969 * t526 - t965 * t529
			+ t965 * t507 - t842 * t529 + t842 * t507 + t910 * t493
			+ t866 * t495 + t872 * t497 - t978 * t481 + t980 * t504
			- t982 * t507 - t984 * t507;
	double t991 = cp2 * t9;
	double t993 = cp2 * b3;
	double t997 = bp2 * b2;
	double t1002 = bp1 * cp3;
	double t1008 = t980 * t511 + t991 * t513 + t993 * t515 - t984 * t517
			- t978 * t517 - t997 * t260 - t997 * t497 + t844 * t481
			- t969 * t194 + t1002 * t101 + t993 * t292 - t872 * t278
			+ 2 * t872 * t276;
	double t1023 = -t945 * t45 - t945 * t66 + t945 * t120 + t945 * t123
			+ t945 * t125 - t945 * t75 - t945 * t51 - t945 * t78 - t1002 * t103
			- t1002 * t106 - t1002 * t109 - t1002 * t111 + t1002 * t60
			+ t1002 * t63;
	double t1039 = t1002 * t45 - t1002 * t118 - t1002 * t120 - t1002 * t123
			- t1002 * t125 + t1002 * t75 + t1002 * t51 + t1002 * t78
			- t945 * t49 - t945 * t57 + t945 * t85 + t945 * t88 - t860 * t111
			+ t860 * t63;
	double t1054 = t860 * t45 - t860 * t118 - t860 * t120 - t860 * t123
			- t860 * t125 + t860 * t75 + t860 * t78 + t1002 * t49 + t1002 * t57
			- t1002 * t85 - t1002 * t88 - t1002 * t91 - t1002 * t33
			+ t1002 * t97;
	double t1063 = ap2 * a3;
	double t1065 = ap2 * a2;
	double t1074 = t1002 * t93 - t905 * t309 - t969 * t191 + t969 * t37
			+ 2 * t969 * t314 - t1063 * t157 - t1065 * t194 - t910 * t196
			- t905 * t199 + 2 * t905 * t201 - t1065 * t162 + t851 * t125
			- t851 * t75;
	double t1089 = -t851 * t51 - t851 * t78 - t851 * t49 - t851 * t57
			+ t851 * t85 + t851 * t88 + t851 * t33 - t851 * t97 - t851 * t93
			- t851 * t101 + t851 * t106 + t851 * t109 - t851 * t63 - t851 * t45;
	double t1108 = -t851 * t66 + t851 * t118 + t851 * t123 + 2 * t905 * t243
			- t997 * t276 - t997 * t278 + t855 * t281 - t855 * t274
			+ 2 * t855 * t286 - t855 * t278 + t877 * t72 - t874 * t25
			+ 2 * t874 * t292 - t870 * t295;
	double t1110 = bp2 * c1;
	double t1119 = bp2 * t18;
	double t1126 = t870 * t297 + 2 * t1110 * t299 - t874 * t72 - t877 * t25
			- t877 * t305 + t1110 * t267 + t851 * t91 - t844 * t3 - t1119 * t7
			- t965 * t11 - t991 * t18 - t842 * t7 - t991 * t20 - t965 * t7;
	double t1136 = cp2 * a1;
	double t1142 = cp2 * c1;
	double t1149 = -t840 * t13 + 2 * t905 * t157 + 2 * t910 * t162
			+ t1136 * t211 - t980 * t211 - t980 * t269 - t872 * t191
			+ t993 * t442 - t1142 * t297 - t1142 * t41 + t1136 * t265
			+ 2 * t1136 * t267 - t1142 * t269;
	double t1165 = -t993 * t451 + 2 * t1136 * t53 - t1136 * t41 - t853 * t57
			+ t853 * t85 + t853 * t88 + t853 * t91 + t853 * t33 - t853 * t97
			- t853 * t93 - t853 * t101 + t853 * t103 - t853 * t60 - t853 * t63;
	double t1173 = ap2 * c1;
	double t1180 = ap2 * a1;
	double t1185 = -t853 * t45 - t853 * t66 + t853 * t123 - t905 * t205
			+ t1063 * t199 - t866 * t211 + 2 * t1173 * t297 - t1173 * t53
			+ 2 * t1173 * t41 + t1173 * t263 - t1180 * t267 - t1173 * t269
			- t1180 * t53 - t910 * t274;
	double t1187 = bp2 * c2;
	double t1204 = -t877 * t183 + 2 * t1187 * t191 + t1187 * t194 - t1187 * t37
			- t1187 * t314 - t874 * t157 - t855 * t194 + 2 * t874 * t309
			- t874 * t179 - t997 * t196 + 2 * t855 * t196 - t855 * t356
			+ t855 * t314 - t877 * t201;
	double t1211 = cp2 * a3;
	double t1222 = t905 * t451 - t910 * t453 - t853 * t49 - t993 * t166
			- t1211 * t170 + 2 * t993 * t172 - t993 * t170 - t993 * t177
			- t993 * t179 - t993 * t183 - t848 * t286 + t848 * t278
			+ t1211 * t25;
	double t1237 = bp2 * b1;
	double t1239 = -t1211 * t292 - t1136 * t297 - t1142 * t299 - t1211 * t303
			- t993 * t25 + 2 * t993 * t305 - t905 * t29 + t945 * t91
			- t1063 * t243 + t997 * t274 + t855 * t577 - t997 * t577
			+ t870 * t580 - t1237 * t504;
	double t1255 = -t1237 * t580 + t846 * t507 + t1119 * t507 - t1237 * t511
			+ t870 * t591 - t842 * t538 + t842 * t517 - t967 * t538
			+ t967 * t517 + t842 * t28 + t967 * t210 + t969 * t554 + t842 * t210
			+ t965 * t481;
	double t1270 = t965 * t28 - t1237 * t591 + t874 * t594 - t877 * t515
			- t877 * t594 + t1119 * t517 + t844 * t517 + t874 * t600
			- t877 * t600 + t846 * t481 + t855 * t260 - t997 * t521
			+ t846 * t161 + t877 * t535;
	double t1287 = t877 * t485 + t877 * t487 + t877 * t540 + t997 * t489
			+ t997 * t546 + t1237 * t491 + t1119 * t161 + t997 * t551
			+ t997 * t493 + t1237 * t495 - t877 * t483 + t910 * t260
			- t1065 * t577;
	double t1302 = t910 * t577 - t1180 * t580 + t866 * t580 - t1180 * t591
			+ t866 * t591 - t1063 * t594 + t905 * t594 - t1063 * t600
			+ t905 * t600 - t1065 * t260 + t905 * t485 + t905 * t487
			+ t910 * t489 + t866 * t491;
	double t1308 = cp2 * c2;
	double t1319 = -t982 * t481 + t872 * t521 + t991 * t523 + t840 * t523
			- t1308 * t526 + t991 * t529 - t991 * t507 + t838 * t529
			- t838 * t507 + t982 * t161 + t993 * t535 + t991 * t28 + t838 * t538
			+ t993 * t540;
	double t1334 = -t838 * t517 + t840 * t538 - t840 * t517 + t838 * t28
			+ t872 * t546 + t984 * t161 + t872 * t551 + t840 * t210
			+ t858 * t101 + t1173 * t272 + t872 * t439 - t838 * t20 - t967 * t11
			- t1119 * t3;
	double t1353 = -t967 * t3 + t1136 * t295 + t858 * t63 + t1002 * t66
			+ t872 * t286 + 2 * t872 * t419 + 2 * t993 * t185 - t993 * t188
			+ t1308 * t554 + t838 * t210 - t991 * t481 + t993 * t483
			+ t965 * t513 - t1308 * t191;
	double t1370 = -t1308 * t314 + 2 * t848 * t194 - t1211 * t179 - t872 * t356
			- t872 * t241 - t872 * t415 - t848 * t356 - t848 * t314
			- t993 * t201 - t872 * t388 + 2 * t993 * t205 - t993 * t426
			+ t980 * t429 - t980 * t299;
	double t1388 = -t993 * t393 - t910 * t276 + 2 * t910 * t278 + t1065 * t274
			- t1065 * t286 - t905 * t72 - t1063 * t292 - t1173 * t299
			+ t1063 * t72 + 2 * t905 * t25 - t905 * t305 - t1173 * t267
			- t855 * t162 - t997 * t388;
	double t1405 = -t877 * t205 - t874 * t199 + t877 * t393 - t870 * t211
			+ 2 * t1237 * t211 - t1110 * t297 + t1110 * t53 - t1110 * t41
			- t1110 * t263 - t870 * t265 - t870 * t267 + 2 * t1110 * t269
			- t1110 * t272 - t870 * t53;
	double t1412 = ap3 * t5;
	double t1415 = ap3 * t1;
	double t1417 = bp3 * t13;
	double t1419 = cp3 * t5;
	double t1421 = bp3 * t18;
	double t1423 = bp1 * ap2;
	double t1425 = bp3 * t20;
	double t1427 = cp3 * a2;
	double t1430 = cp3 * t9;
	double t1433 = cp3 * t1;
	double t1435 = -t1412 * t3 - t1412 * t7 - t1415 * t3 - t1417 * t7
			- t1419 * t13 - t1421 * t7 - t1423 * t88 - t1425 * t3 - t1427 * t281
			- t1415 * t11 - t1430 * t18 - t1419 * t20 - t1433 * t13;
	double t1436 = bp2 * cp1;
	double t1439 = ap3 * b1;
	double t1443 = ap2 * cp1;
	double t1445 = bp3 * b3;
	double t1447 = bp3 * a3;
	double t1450 = cp3 * a1;
	double t1452 = cp3 * b2;
	double t1455 = ap3 * c1;
	double t1459 = bp1 * cp2;
	double t1461 = -t1436 * t88 + t1423 * t45 + t1439 * t299 + t1436 * t60
			+ t1436 * t51 + t1443 * t33 + t1445 * t72 - t1447 * t170
			- t1425 * t11 + t1450 * t295 + 2 * t1452 * t276 + 2 * t1455 * t41
			+ t1450 * t211 + t1459 * t120;
	double t1463 = bp3 * a1;
	double t1465 = ap1 * cp2;
	double t1469 = bp2 * ap1;
	double t1477 = bp3 * a2;
	double t1481 = -t1463 * t53 + t1465 * t57 + t1423 * t51 + t1436 * t45
			+ t1469 * t123 + t1443 * t88 - t1459 * t97 - t1423 * t118
			+ t1469 * t120 + t1469 * t33 - t1452 * t191 + t1477 * t314
			+ t1443 * t85 + t1459 * t111;
	double t1484 = cp3 * a3;
	double t1488 = cp3 * b3;
	double t1492 = cp3 * b1;
	double t1498 = ap3 * b3;
	double t1500 = -t1436 * t33 + t1459 * t123 + t1484 * t25 - t1459 * t57
			- t1477 * t356 + t1488 * t442 - t1488 * t451 + t1452 * t453
			- t1492 * t269 + t1488 * t292 - t1488 * t426 + t1452 * t286
			- t1452 * t415 - t1498 * t179;
	double t1504 = ap3 * b2;
	double t1517 = -t1498 * t188 - t1504 * t356 - t1504 * t241 - t1504 * t439
			- t1498 * t442 + t1498 * t451 - t1504 * t453 + t1439 * t269
			- t1498 * t292 - t1439 * t429 - t1459 * t63 - t1459 * t49
			- t1459 * t60;
	double t1532 = -t1459 * t66 + t1459 * t33 - t1459 * t101 + t1459 * t85
			+ t1459 * t109 - t1459 * t51 + t1459 * t103 + t1459 * t106
			- t1459 * t75 + t1459 * t91 - t1459 * t78 + t1423 * t97
			+ t1423 * t63 + t1423 * t49;
	double t1548 = t1423 * t60 + t1423 * t66 - t1423 * t33 + t1423 * t101
			- t1423 * t123 - t1423 * t125 + t1423 * t93 - t1423 * t85
			- t1423 * t111 - t1423 * t109 - t1423 * t120 + t1423 * t57
			- t1423 * t103 - t1423 * t106;
	double t1563 = t1423 * t75 - t1423 * t91 + t1423 * t78 + t1465 * t49
			+ t1465 * t101 - t1465 * t123 - t1465 * t125 + t1465 * t93
			- t1465 * t85 - t1465 * t111 - t1465 * t109 - t1465 * t120
			- t1465 * t118 + t1465 * t51;
	double t1580 = t1465 * t45 - t1465 * t103 - t1465 * t106 + t1465 * t75
			+ t1465 * t78 + t1436 * t97 + t1436 * t63 + t1436 * t49
			+ t1436 * t66 + t1436 * t101 - t1436 * t123 - t1436 * t125
			- t1436 * t111;
	double t1595 = -t1436 * t109 - t1436 * t120 - t1436 * t118 + t1436 * t57
			- t1436 * t103 + t1450 * t265 + t1459 * t88 + t1459 * t118
			- t1465 * t88 + t1504 * t415 + t1459 * t125 + t1465 * t97
			+ t1465 * t63 - t1477 * t162;
	double t1611 = -t1436 * t106 + t1436 * t75 - t1436 * t91 + t1436 * t78
			- t1469 * t97 - t1469 * t63 - t1469 * t49 - t1469 * t60
			- t1469 * t101 + t1469 * t125 - t1469 * t93 + t1469 * t85
			+ t1469 * t111 + t1469 * t109;
	double t1626 = t1469 * t88 + t1469 * t118 - t1469 * t51 - t1469 * t57
			- t1469 * t45 - t1469 * t75 + t1469 * t91 - t1469 * t78
			- t1443 * t97 - t1443 * t63 - t1443 * t49 - t1443 * t60
			- t1443 * t66 - t1443 * t101;
	double t1643 = t1443 * t123 + t1443 * t125 + t1443 * t118 - t1443 * t51
			- t1443 * t57 - t1443 * t45 + t1443 * t103 + t1443 * t106
			- t1443 * t75 + t1443 * t91 - t1443 * t78 + 2 * t1504 * t162
			- t1484 * t170;
	double t1649 = ap3 * a3;
	double t1661 = -t1488 * t188 - t1445 * t29 - t1498 * t29 - t1452 * t278
			- t1427 * t314 + t1649 * t72 - t1488 * t166 + 2 * t1488 * t172
			- t1504 * t286 + 2 * t1427 * t194 - t1498 * t170 - t1445 * t183
			- t1445 * t309 - t1443 * t93;
	double t1669 = bp3 * b2;
	double t1674 = bp3 * c2;
	double t1682 = t1443 * t111 + 2 * t1447 * t29 - t1447 * t243 - t1445 * t177
			- t1447 * t179 - t1669 * t196 + 2 * t1477 * t196 - t1669 * t388
			+ 2 * t1674 * t191 + t1674 * t194 - t1674 * t314 - t1674 * t37
			- t1477 * t346 - t1477 * t194;
	double t1687 = bp3 * b1;
	double t1691 = bp3 * c1;
	double t1701 = t1445 * t199 - t1445 * t205 - t1445 * t201 - t1447 * t199
			+ 2 * t1687 * t211 - t1463 * t211 + t1691 * t53 - t1691 * t297
			- t1691 * t272 + 2 * t1691 * t269 - t1463 * t295 - t1447 * t25
			- t1445 * t305 + t1447 * t303;
	double t1721 = -t1447 * t72 + 2 * t1691 * t299 + t1691 * t267 - t1691 * t41
			- t1691 * t263 - t1463 * t265 - t1463 * t267 + t1463 * t41
			+ t1669 * t274 - t1669 * t276 + t1477 * t281 - t1477 * t274
			+ 2 * t1477 * t286;
	double t1725 = ap3 * t9;
	double t1734 = cp3 * c2;
	double t1738 = -t1477 * t278 - t1649 * t243 - t1421 * t3 - t1725 * t7
			- t1430 * t20 - t1725 * t11 - t1436 * t85 + t1465 * t60
			+ t1465 * t66 - t1459 * t45 + t1504 * t191 - t1734 * t191
			+ t1443 * t109 + t1443 * t120;
	double t1750 = cp3 * t13;
	double t1752 = cp3 * t20;
	double t1756 = -t1465 * t33 - t1459 * t93 - t1465 * t91 + t1436 * t93
			+ t1427 * t346 + t1463 * t297 + t1430 * t513 - t1419 * t507
			+ t1419 * t529 + t1430 * t523 - t1750 * t507 - t1752 * t481
			+ t1750 * t161 + t1433 * t523;
	double t1768 = cp3 * t18;
	double t1772 = t1488 * t535 + t1492 * t504 + t1430 * t28 - t1419 * t517
			+ t1488 * t540 + t1492 * t511 + t1433 * t538 + t1419 * t28
			- t1433 * t517 + t1452 * t546 + t1452 * t521 + t1768 * t161
			- t1734 * t526 + t1452 * t551;
	double t1788 = t1430 * t529 + t1433 * t210 + t1734 * t554 + t1419 * t538
			+ t1419 * t210 + t1452 * t497 + t1488 * t515 - t1750 * t481
			- t1752 * t517 - t1768 * t517 + t1725 * t507 + t1725 * t481
			+ t1725 * t513;
	double t1798 = ap3 * c2;
	double t1804 = t1412 * t507 - t1412 * t529 - t1725 * t523 - t1415 * t523
			+ t1725 * t28 + t1412 * t517 - t1415 * t538 + t1412 * t28
			+ t1415 * t517 + t1798 * t526 - t1725 * t529 + t1415 * t210
			+ t1798 * t554 - t1412 * t538;
	double t1811 = ap3 * a2;
	double t1824 = t1412 * t210 - t1498 * t309 - t1649 * t157 + 2 * t1504 * t278
			- t1811 * t162 - t1798 * t191 - t1798 * t194 + 2 * t1798 * t314
			+ t1798 * t37 - t1811 * t194 - t1498 * t205 + 2 * t1498 * t201
			+ t1649 * t199 - t1439 * t211;
	double t1830 = ap3 * a1;
	double t1842 = -t1455 * t53 + 2 * t1455 * t297 + t1455 * t272 - t1455 * t269
			- t1830 * t53 - t1649 * t292 - t1498 * t72 - t1498 * t305
			+ 2 * t1498 * t25 - t1455 * t299 - t1455 * t267 + t1455 * t263
			- t1830 * t267 - t1504 * t274;
	double t1860 = -t1504 * t276 + t1811 * t274 - t1811 * t286 - t1488 * t170
			- t1488 * t177 - t1488 * t179 - t1488 * t183 + 2 * t1488 * t185
			- t1484 * t179 - t1452 * t356 - t1452 * t241 - t1427 * t356
			- t1452 * t388;
	double t1865 = cp3 * c1;
	double t1878 = -t1734 * t314 - t1488 * t201 - t1488 * t393 - t1492 * t211
			- t1865 * t297 - t1865 * t269 + 2 * t1450 * t53 - t1450 * t297
			- t1484 * t292 + 2 * t1488 * t305 - t1488 * t25 - t1484 * t303
			- t1865 * t299 - t1865 * t41;
	double t1895 = 2 * t1450 * t267 - t1450 * t41 + t1452 * t439 - t1417 * t11
			- t1433 * t18 + t1447 * t594 + t1463 * t591 + t1421 * t507
			- t1445 * t483 - t1669 * t260 - t1669 * t577 - t1687 * t580
			+ t1417 * t507 + t1425 * t481;
	double t1910 = t1417 * t161 + t1445 * t535 - t1687 * t504 + t1445 * t485
			- t1687 * t591 + t1447 * t600 + t1445 * t487 + t1445 * t540
			- t1687 * t511 + t1477 * t577 + t1669 * t489 + t1669 * t546
			+ t1687 * t491 - t1445 * t594;
	double t1927 = -t1669 * t521 + t1421 * t161 + t1669 * t551 + t1669 * t493
			+ t1463 * t580 - t1669 * t497 - t1445 * t600 - t1445 * t515
			+ t1477 * t260 + t1687 * t495 + t1417 * t481 + t1425 * t517
			+ t1421 * t517 - t1649 * t594;
	double t1942 = -t1830 * t591 + t1504 * t260 + t1504 * t577 + t1439 * t580
			+ t1498 * t485 + t1439 * t591 - t1649 * t600 + t1498 * t487
			- t1811 * t577 + t1504 * t489 + t1439 * t491 + t1498 * t594
			+ t1504 * t493 - t1830 * t580;
	double t1961 = t1498 * t600 - t1811 * t260 + t1439 * t495 - t1430 * t507
			- t1430 * t481 - t1768 * t507 + t1488 * t483 - t1469 * t66
			+ t1492 * t429 + t1427 * t278 + 2 * t1447 * t309 + 2 * t1452 * t419
			- t1447 * t157 + 2 * t1447 * t292;
	double t1979 = t1445 * t393 - t1498 * t166 - t1498 * t199 - t1427 * t286
			- t1445 * t25 + t1498 * t426 - t1669 * t278 - t1492 * t299
			+ t1469 * t103 + t1469 * t106 + 2 * t1498 * t157 - t1504 * t196
			+ 2 * t1488 * t205 + 2 * t1498 * t243;
	cg3[0] = -(t698 + t713 + t667 + t47 + t82 + t150 + t181 + t114 + t134 + t215
			+ t230 + t253 + t283 + t311 + t328 + t348 + t370 + t390 + t410
			+ t436 + t458 + t476 + t501 + t528 + t548 + t566 + t583 + t603
			+ t618 + t634 + t649 + t682) * t836;
	cg3[1] = -(t1023 + t1185 + t1270 + t862 + t1126 + t1353 + t1255 + t1074
			+ t1287 + t1149 + t882 + t986 + t1319 + t1388 + t1089 + t1204 + t932
			+ t1405 + t1222 + t898 + t1165 + t1302 + t1334 + t1008 + t1370
			+ t1039 + t964 + t1108 + t1054 + t948 + t1239 + t916) * t836;
	cg3[2] = -(t1548 + t1788 + t1942 + t1626 + t1772 + t1961 + t1643 + t1860
			+ t1804 + t1842 + t1661 + t1701 + t1721 + t1895 + t1481 + t1435
			+ t1611 + t1682 + t1595 + t1563 + t1580 + t1532 + t1517 + t1738
			+ t1461 + t1500 + t1910 + t1756 + t1878 + t1979 + t1824 + t1927)
			* t836;

}

/**
 *	2---3    +-1-+   /\      1b--2b
 *	|   |    2   3  4  5      | X |
 *	0---1    +-0-+ /    \    1a--2a
 *	Input l[i][j] and d[i][j]=l[i][j]^2
 */
void Utils::lenToTetra(double* lengths, double* loc1a, double* loc2a,
		double* loc1b, double* loc2b, bool positive) {

//	Output (four points) q
	double l[3][4];
	double d[3][4];
	l[0][1] = lengths[0];
	l[2][3] = lengths[1];
	l[0][2] = lengths[2];
	l[1][3] = lengths[3];
	l[0][3] = lengths[4];
	l[1][2] = lengths[5];

	d[0][1] = l[0][1] * l[0][1];
	d[2][3] = l[2][3] * l[2][3];
	d[0][2] = l[0][2] * l[0][2];
	d[1][3] = l[1][3] * l[1][3];
	d[0][3] = l[0][3] * l[0][3];
	d[1][2] = l[1][2] * l[1][2];

	double t2 = 1.0 / l[0][1];
	double t5 = pow(d[0][1], 2);
	double t6 = d[0][2] * d[0][1];
	double t8 = d[0][1] * d[1][2];
	double t10 = pow(d[0][2], 2);
	double t11 = d[0][2] * d[1][2];
	double t13 = pow(d[1][2], 2);
	double t14 = t5 - 2 * t6 - 2 * t8 + t10 - 2 * t11 + t13;
	double t15 = pow(l[0][1], 2);
	double t18 = sqrt(-1.0 / t15 * t14);
	double t23 = d[1][3] * d[0][1];
	double t24 = d[0][3] * d[0][1];
	double t25 = d[0][3] * d[0][2];
	double t26 = d[0][3] * d[1][2];
	double t27 = d[1][3] * d[0][2];
	double t28 = d[1][3] * d[1][2];
	double t29 = d[2][3] * d[0][1];
	double t33 = sqrt(-d[0][1] * t14);
	double t39 = pow(d[1][3], 2);
	double t42 = pow(d[0][3], 2);
	double t50 = d[1][3] * t24 - d[0][2] * t23 + t39 * d[0][2] - d[1][2] * t24
			+ t42 * d[1][2] + d[1][3] * t10 + d[1][2] * t6 + d[0][3] * t13
			- d[1][3] * t25 - d[1][3] * t26 - d[0][3] * t11;
	double t61 = pow(d[2][3], 2);
	double t63 = -d[1][3] * t11 + d[2][3] * t25 - d[2][3] * t26 - d[2][3] * t27
			+ d[2][3] * t28 - d[2][3] * t23 - d[2][3] * t24 - d[0][2] * t29
			- d[1][2] * t29 + d[2][3] * t5 + t61 * d[0][1];
	double t67 = sqrt(1.0 / t14 * (t50 + t63));

	//output 0,1,2 = x,y,z
//	 * helix a dumbbell = q0 and q1
	loc1a[0] = 0;
	loc1a[1] = 0;
	loc1a[2] = 0;
	loc2a[0] = l[0][1];
	loc2a[1] = 0;
	loc2a[2] = 0;
//	 * helix b dumbbell = q2 and q3
	loc1b[0] = t2 * (d[0][1] + d[0][2] - d[1][2]) / 2;
	loc1b[1] = t18 / 2;
	loc1b[2] = 0;
	loc2b[0] = t2 * (d[0][1] + d[0][3] - d[1][3]) / 2;
	loc2b[1] = 1.0 / t33
			* (-t5 + t23 + t24 - t25 + t26 + t27 - t28 - 2 * t29 + t6 + t8) / 2;
	if (positive) {	//plus or minus
		loc2b[2] = t67;
	} else {
		loc2b[2] = -1 * t67;
	}
}

vector<double> Utils::getIdentityMatrix() {
	static double eye[] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 };
	return vector<double>(eye, eye + 12);
}


vector<double> Utils::getTransMatrix(double from[][3], double to[][3]) {
	vector<double> matrix;
	double o[] = { 0, 0, 0 };
	double x[] = { 1, 0, 0 };
	double y[] = { 0, 1, 0 };
	double z[] = { 0, 0, 1 };
	double no[3], nx[3], ny[3], nz[3];

	Utils::matApp(o, from[0], from[1], from[2], to[0], to[1], to[2], no);
	Utils::matApp(x, from[0], from[1], from[2], to[0], to[1], to[2], nx);
	Utils::matApp(y, from[0], from[1], from[2], to[0], to[1], to[2], ny);
	Utils::matApp(z, from[0], from[1], from[2], to[0], to[1], to[2], nz);
	// row major matrix
	for (int i = 0; i < 3; i++) {
		matrix.push_back(nx[i] - no[i]);
		matrix.push_back(ny[i] - no[i]);
		matrix.push_back(nz[i] - no[i]);
		matrix.push_back(no[i]);
	}

	return matrix;

}

vector<double> Utils::getSingleTransformFromTwo(const vector<double> &ta,
		const vector<double> &tb) {
	Matrix3d ra, rb;	// rotation matrix
	Vector3d tla, tlb; // translation vector

	ra << ta[0], ta[1], ta[2], ta[4], ta[5], ta[6], ta[8], ta[9], ta[10];
	rb << tb[0], tb[1], tb[2], tb[4], tb[5], tb[6], tb[8], tb[9], tb[10];

	tla << ta[3], ta[7], ta[11];
	tlb << tb[3], tb[7], tb[11];

	// fix a instead of b
	Matrix3d newr = ra.inverse() * rb;
	Vector3d newt = ra.inverse() * (tlb - tla);


	vector<double> res;
	for (int i = 0; i < 3; i++) {
		res.push_back(newr(i, 0));
		res.push_back(newr(i, 1));
		res.push_back(newr(i, 2));
		res.push_back(newt(i));
	}
	return res;
}

vector<double> Utils::multiTrans(const vector<double> &ta,
		const vector<double> &tb) {
	Matrix3d ra, rb;	// rotation matrix
	Vector3d tla, tlb; // translation vector

	ra << ta[0], ta[1], ta[2], ta[4], ta[5], ta[6], ta[8], ta[9], ta[10];
	rb << tb[0], tb[1], tb[2], tb[4], tb[5], tb[6], tb[8], tb[9], tb[10];

	tla << ta[3], ta[7], ta[11];
	tlb << tb[3], tb[7], tb[11];

	// apply a and then b
	Matrix3d newr = rb * ra;
	Vector3d newt = rb * tla + tlb;

	vector<double> res;
	for (int i = 0; i < 3; i++) {
		res.push_back(newr(i, 0));
		res.push_back(newr(i, 1));
		res.push_back(newr(i, 2));
		res.push_back(newt(i));
	}
	return res;
}

vector<double> Utils::getMatrixFromFileFromTo(const char *filename) {
	double true_from[3][3], true_to[3][3];
	ifstream true_file(filename);
	if (true_file.fail()) {
		cout << "read fromto fail" << endl;
		true_file.close();
		return Utils::getIdentityMatrix();
	}
	for (int i = 0; i < 3; i++) {
		true_file >> true_from[i][0] >> true_from[i][1] >> true_from[i][2];
		true_file >> true_to[i][0] >> true_to[i][1] >> true_to[i][2];
	}
	true_file.close();
	return Utils::getTransMatrix(true_from, true_to);

}

vector<vector<double> > Utils::getMatrixFromFileMat(const char* filename) {
	vector<vector<double> > res;
	ifstream true_file(filename);
	if (true_file.fail()) {
		cout << "read matfile fail" << endl;
		return res;
	}
	while (true_file.good()) {
		vector<double> mat(12);
		for (int i = 0; i < 12; i++) {
			if (!true_file.good() || true_file.eof()) {
				true_file.close();
				return res;
			}

			true_file >> mat[i];
		}
		res.push_back(mat);
	}
	true_file.close();
	return res;

}
