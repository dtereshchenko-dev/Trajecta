#pragma once

#include <QSignalSpy>
#include <QTest>

#include "geometry/OrbitGeometry.h"

#include <cmath>
#include <cstring>
#include <limits>

class TestCometOrbits : public QObject {
	Q_OBJECT

	static void verifyFiniteGeometry(const OrbitGeometry& geometry)
	{
		const QByteArray vertices = geometry.vertexData();
		QVERIFY(!vertices.isEmpty());
		QCOMPARE(vertices.size() % qsizetype(3 * sizeof(float)), qsizetype{0});
		const QVector3D minimum = geometry.boundsMin();
		const QVector3D maximum = geometry.boundsMax();
		for (int axis = 0; axis < 3; ++axis) {
			QVERIFY(std::isfinite(minimum[axis]));
			QVERIFY(std::isfinite(maximum[axis]));
			QVERIFY(minimum[axis] <= maximum[axis]);
		}
		for (qsizetype offset = 0; offset < vertices.size(); offset += 3 * sizeof(float)) {
			float point[3];
			std::memcpy(point, vertices.constData() + offset, sizeof(point));
			for (int axis = 0; axis < 3; ++axis) {
				QVERIFY(std::isfinite(point[axis]));
				QVERIFY(point[axis] >= minimum[axis]);
				QVERIFY(point[axis] <= maximum[axis]);
			}
			const double radius = std::hypot(static_cast<double>(point[0]), static_cast<double>(point[1]),
			                                 static_cast<double>(point[2]));
			QVERIFY(radius <= geometry.maximumRadius() * (1.0 + 1e-6));
		}
	}

   private slots:

	void perigeePropertyDefaultAndNotifications()
	{
		OrbitGeometry geometry;
		QCOMPARE(geometry.perigeeAu(), 0.0);
		QSignalSpy spy(&geometry, &OrbitGeometry::geometryChanged);
		geometry.setPerigeeAu(10.9);
		QCOMPARE(geometry.perigeeAu(), 10.9);
		QCOMPARE(spy.count(), 1);
		geometry.setPerigeeAu(10.9);
		QCOMPARE(spy.count(), 1);
	}

	void finiteGeometryAcrossParabolicBoundary()
	{
		for (const double e : {0.998, std::nextafter(1.0, 0.0), 1.0, std::nextafter(1.0, 2.0), 1.2, 3.3769}) {
			OrbitGeometry geometry;
			geometry.setPerigeeAu(10.9);
			geometry.setEccentricity(e);
			geometry.setSemiMajorAxis(e == 1.0 ? std::numeric_limits<double>::infinity() : 10.9 / (1.0 - e));
			geometry.setInclination(95.538);
			geometry.setAscendingNode(190.066);
			geometry.setArgPericentre(326.375);
			verifyFiniteGeometry(geometry);
			QVERIFY(std::abs(geometry.perihelionPos().length() - 10.9f) < 1.0e-4f);
			const int vertexCount = geometry.vertexData().size() / (3 * sizeof(float));
			QCOMPARE(vertexCount, e < 1.0 ? geometry.segments() + 1 : geometry.segments());
		}
	}

	void qGeometryDoesNotDependOnSemimajorAxis()
	{
		OrbitGeometry geometry;
		geometry.setPerigeeAu(2.0);
		geometry.setEccentricity(1.0);
		const auto expected = geometry.vertexData();
		for (const double a :
		     {0.0, 1.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
			geometry.setSemiMajorAxis(a);
			QCOMPARE(geometry.vertexData(), expected);
			verifyFiniteGeometry(geometry);
		}
	}

	void highEccentricityExtentCoversAphelion()
	{
		OrbitGeometry geometry;
		geometry.setSegments(512);
		geometry.setPerigeeAu(2.2666);
		geometry.setEccentricity(0.993748556);
		const double aphelion = geometry.perigeeAu() * (1.0 + geometry.eccentricity()) /
		                        (1.0 - geometry.eccentricity());
		QVERIFY(std::abs(geometry.maximumRadius() - aphelion) < 1e-8 * aphelion);
		verifyFiniteGeometry(geometry);
	}

	void highEccentricityEllipseHasSmoothEnds()
	{
		OrbitGeometry geometry;
		geometry.setSegments(512);
		geometry.setInclination(0.0);
		geometry.setAscendingNode(0.0);
		geometry.setArgPericentre(0.0);
		const double a = 345.003531;
		const double e = 0.999146;
		geometry.setSemiMajorAxis(a);
		geometry.setEccentricity(e);
		geometry.setPerigeeAu(a * (1.0 - e));
		const QByteArray vertices = geometry.vertexData();
		const auto point = [&](int index) {
			float coords[3];
			std::memcpy(coords, vertices.constData() + index * 3 * sizeof(float), sizeof(coords));
			return QVector3D(coords[0], coords[1], coords[2]);
		};
		QCOMPARE(point(0), point(512));
		const auto smoothTurn = [](const QVector3D& before, const QVector3D& at, const QVector3D& after) {
			const QVector3D incoming = (at - before).normalized();
			const QVector3D outgoing = (after - at).normalized();
			return QVector3D::dotProduct(incoming, outgoing);
		};
		QVERIFY(smoothTurn(point(255), point(256), point(257)) > 0.99f);
		QVERIFY(smoothTurn(point(511), point(0), point(1)) > 0.99f);
		const double b = a * std::sqrt((1.0 - e) * (1.0 + e));
		for (int i = 0; i <= 512; ++i) {
			const auto vertex = point(i);
			const double u = (vertex.x() + a * e) / a;
			const double v = vertex.z() / b;
			QVERIFY(std::abs(u * u + v * v - 1.0) < 1e-5);
		}
	}

	void missingPerigeeKeepsConventionalGeometry()
	{
		OrbitGeometry conventional;
		conventional.setEccentricity(0.5);
		conventional.setSemiMajorAxis(2.0);
		OrbitGeometry withPerigee;
		withPerigee.setEccentricity(0.5);
		withPerigee.setSemiMajorAxis(2.0);
		withPerigee.setPerigeeAu(1.0);
		QCOMPARE(conventional.vertexData(), withPerigee.vertexData());
		QCOMPARE(conventional.indexData(), withPerigee.indexData());
		verifyFiniteGeometry(conventional);
	}

	void invalidGeometryClearsVerticesAndPerihelion()
	{
		OrbitGeometry geometry;
		geometry.setPerigeeAu(2.0);
		geometry.setEccentricity(1.0);
		geometry.setSemiMajorAxis(std::numeric_limits<double>::infinity());
		QVERIFY(!geometry.vertexData().isEmpty());
		geometry.setPerigeeAu(0.0);
		QVERIFY(geometry.vertexData().isEmpty());
		QCOMPARE(geometry.maximumRadius(), 0.0);
		QCOMPARE(geometry.perihelionPos(), QVector3D());
		for (const double q :
		     {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
			geometry.setPerigeeAu(q);
			QVERIFY(geometry.vertexData().isEmpty());
		}
		geometry.setPerigeeAu(2.0);
		geometry.setEccentricity(-0.1);
		QVERIFY(geometry.vertexData().isEmpty());
		geometry.setEccentricity(1.0);
		geometry.setSegments(0);
		QVERIFY(geometry.vertexData().isEmpty());
		QCOMPARE(geometry.perihelionPos(), QVector3D());
	}
};
