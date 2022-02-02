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
#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

using Eigen::Vector3d;
using Eigen::Matrix3d;

using namespace std;

class Vector3 {
public:
	Vector3() {
		Vector3(0, 0, 0);
	}

	Vector3(double x, double y, double z) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}

	Vector3(double *p) {
		coord[0] = p[0];
		coord[1] = p[1];
		coord[2] = p[2];
	}

	Vector3 operator+(const Vector3& v) {
		return Vector3(coord[0] + v.coord[0], coord[1] + v.coord[1],
				coord[2] + v.coord[2]);
	}

	Vector3 operator-(const Vector3& v) {
		return Vector3(coord[0] - v.coord[0], coord[1] - v.coord[1],
				coord[2] - v.coord[2]);
	}

	Vector3 operator*(double k) {
		return Vector3(coord[0] * k, coord[1] * k, coord[2] * k);
	}

	Vector3 operator/(double k) {
		return Vector3(coord[0] / k, coord[1] / k, coord[2] / k);
	}

	Vector3 normalize() {
		double len = sqrt(this->squared_length());
		return (*this) / len;
	}
	double squared_length() {
		return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
	}

	double x() {
		return coord[0];
	}
	double y() {
		return coord[1];
	}
	double z() {
		return coord[2];
	}

	double& operator[](int idx) {
		return coord[idx];
	}

	double* getLoc() {
		return coord;
	}

	Vector3 trans(const std::vector<double>& t) {
		double ncoord[3];
		for (int i = 0; i < 3; i++) {
			ncoord[i] = 0;

			// rotation
			for (int j = 0; j < 3; j++) {
				ncoord[i] += t[i * 4 + j] * coord[j];
			}

			// translation
			ncoord[i] += t[i * 4 + 3];
		}
		return Vector3(ncoord);
	}

private:
	double coord[3];
};

typedef Vector3 Vector;

/*
 * Some functions to help out with 3D problems.
 */
class Utils {
public:

	/** @return The index of mypair in myvector if found, -1 otherwise */
	static int findPair(vector<pair<int, int> > myvector,
			pair<int, int> mypair);

	/**
	 * @return The distance in L2 norm between a and b
	 * @param a b 3 dimensional arrays
	 */
	static double dist(double *a, double *b);
	static double dist(Vector3d a, Vector3d b);

	/**
	 * @return The distance in L2 norm between a and b
	 * @param a b 6 dimensional arrays
	 */
	static double dist6(double *a, double *b);

	/**
	 * @return The magnitude (norm) of input a
	 * @param a 3 dimensional array
	 */
	static double mag(double *a);

	/*
	 * @param hue A value [0..1]
	 * @param[out] color The 3 part RGB array
	 */
	static void huetoColor(double hue, float color[3]);

	/**
	 * @return 3 dimensional array that is the cross product between b and c
	 * @param b c 3 dimensional arrays
	 */
	static double* crossProd(double *b, double *c);

	/**
	 * @return A scalar value that is the dot product between a and b
	 * @param a b 3 dimensional arrays
	 */
	static double dotProd(double *a, double *b);

	/**
	 * @return 3 dimensional array that is the scalar product between b and c
	 * @param b Scalar
	 * @param c 3 dimensional array
	 */
	static double* scalProd(double b, double* c);

	/**
	 * @return 3 dimensional array that is the b minus c
	 * @param b c 3 dimensional arrays
	 */
	static double* vectSub(double *b, double *c);

	/**
	 3---4    +-1-+   /\
			|   |    2   3  4  5
	 1---2    +-0-+ /    \
	 *
	 * @brief Computes the volume of tetrahedron
	 *
	 * This method sometimes give negative value (may be huge neg value)  when volume is close to zero. It cannot handle critic cases
	 *
	 * @param len 6dimensional array that holds the lengths of edges of the tetrahedron
	 */
	static double volumeTetra(double *len);

	/**
	 * @brief Computes the volume of tetrahedron (Heron-type formula)
	 * This function always give nan value for negative values, because of sqrt.
	 * @param len 6dimensional array that holds the lengths of edges of the tetrahedron
	 */
	static double volumeTetra2(double *len);

	/*
	 * @param x is a 3 part array describing the point to be moved.
	 * @param a b c  are the original locations
	 * @param ap bp cp are the destination locations
	 * @param[out] cg3 The transformed location of x
	 */
	static void matApp(double* x, double* a, double* b, double* c, double* ap,
			double* bp, double* cp, double cg3[3]);

	/*
	 * @param x is a 3 part array describing the point to be moved.
	 * @param a b c  are the original locations
	 * @param ap bp cp are the destination locations
	 * @return The transformed location of x
	 */
	static double* matApp2(double* x, double* a, double* b, double* c,
			double* ap, double* bp, double* cp);

	/*
	 * @param x is a 3 part array describing the point to be moved.
	 * @param a b c  are the original locations
	 * @param ap bp cp are the destination locations
	 * @return The transformed location of x
	 */
	static double* newMatApp(double* x, double* a, double* b, double* c,
			double* ap, double* bp, double* cp);

	/**
	 * @brief Computes the location of a point that lies in the same plane with a b and c
	 * @param a b c  are the locations
	 * @param ap bp cp are the lengths of the edges from unknown point x to those 3 points a b c
	 * @return 3 dimensional array
	 */
	static double* find4thPoint_inTheSamePlane(double* a, double* b, double* c,
			double ap, double bp, double cp);

	/**
	 * @brief This method takes a set of 6 lengths and outputs the locations of the 4 points
	 *
	 * @param lengths 6dimensional array that holds the lengths of edges of the tetrahedron
	 * @param positive Determines which side of the z-axis to put 4th point.
	 *
	 * @param[out] loc1a First point located at origin
	 * @param[out] loc2a Second point located on x axis
	 * @param[out] loc1b Third point located on x-y plane
	 * @param[out] loc2b Fourth point located on x-y-z plane
	 */
	static void lenToTetra(double* lengths, double* loc1a, double* loc2a,
			double* loc1b, double* loc2b, bool positive);

	/**
	 * @brief This method returns a matrix that is combination of rotation and translation.
	 *
	 * @param fa Original Cartesian coordinates of three points
	 * @param ta Destination Cartesian coordinates of three points
	 *
	 * todo For better perfomance try to increase speed of this function
	 */
	static std::vector<double> getTransMatrix(double fa[][3], double ta[][3]);

	/**
	 * @param ta Transformation matrix that is combination of rotation and translation.
	 * @param tb Transformation matrix that is combination of rotation and translation.
	 * @return Transformation matrix that is combination of ta and tb
	 */
	static std::vector<double> getSingleTransformFromTwo(
			const std::vector<double> &ta, const std::vector<double> &tb);

	/**
	 * @param ta Transformation matrix that is combination of rotation and translation.
	 * @param tb Transformation matrix that is combination of rotation and translation.
	 * @return Transformation matrix that is combination of ta and tb (First apply a and then b)
	 */
	static std::vector<double> multiTrans(const std::vector<double> &ta,
			const std::vector<double> &tb);

	/** @return Identity matrix */
	static std::vector<double> getIdentityMatrix();

	/**
	 * @brief This method returns a matrix that is combination of rotation and translation.
	 *
	 * @param filename The path of the file that holds original Cartesian coordinates of three points
	 * and destination Cartesian coordinates of three points
	 *
	 * @see getTransMatrix
	 */
	static std::vector<double> getMatrixFromFileFromTo(const char* filename);

	/**
	 * @brief This method reads and returns the matrix with 12 elements.
	 *
	 * @param filename The path of the file that holds the matrix
	 */
	static std::vector<std::vector<double> > getMatrixFromFileMat(
			const char* filename);

	/** @return The sign of the number no */
	static int sign(double no);

	/** @return The absolute value of the number no */
	static double abs(double no);

	/** @return The acos value of the number no */
	static double acoss(double no);

	/**
	 * @brief COnverts rotation matrix to eular angles
	 * @param rmat Rotation matrix
	 * @return A vector composed of (phi, cos_theta, psi)
	 */
	static Vector3d RotMatToEuler(Matrix3d rmat);
};

#endif /* UTILS_H_ */
