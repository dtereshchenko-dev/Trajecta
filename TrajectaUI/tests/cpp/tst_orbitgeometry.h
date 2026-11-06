#ifndef TST_ORBITGEOMETRY_H
#define TST_ORBITGEOMETRY_H

#include <QTest>
#include <QSignalSpy>
#include <QVector3D>
#include <cmath>
#include "geometry/OrbitGeometry.h"

class TestOrbitGeometry : public QObject {
	Q_OBJECT

   private slots:

	// --- Default values ---

	void defaultValues()
	{
		OrbitGeometry geom;
		QCOMPARE(geom.eccentricity(), 3.3769);
		QCOMPARE(geom.semiMajorAxis(), -0.8449);
		QCOMPARE(geom.segments(), 200);
	}

	// --- Setters emit signal ---

	void setEccentricity_emitsSignal()
	{
		OrbitGeometry geom;
		QSignalSpy spy(&geom, &OrbitGeometry::geometryChanged);
		geom.setEccentricity(0.5);
		QCOMPARE(geom.eccentricity(), 0.5);
		QVERIFY(spy.count() >= 1);
	}

	void setSemiMajorAxis_emitsSignal()
	{
		OrbitGeometry geom;
		QSignalSpy spy(&geom, &OrbitGeometry::geometryChanged);
		geom.setSemiMajorAxis(2.0);
		QCOMPARE(geom.semiMajorAxis(), 2.0);
		QVERIFY(spy.count() >= 1);
	}

	void setInclination_emitsSignal()
	{
		OrbitGeometry geom;
		QSignalSpy spy(&geom, &OrbitGeometry::geometryChanged);
		geom.setInclination(30.0);
		QCOMPARE(geom.inclination(), 30.0);
		QVERIFY(spy.count() >= 1);
	}

	void setSegments_emitsSignal()
	{
		OrbitGeometry geom;
		QSignalSpy spy(&geom, &OrbitGeometry::geometryChanged);
		geom.setSegments(100);
		QCOMPARE(geom.segments(), 100);
		QVERIFY(spy.count() >= 1);
	}

	// --- No signal on same value ---

	void setSameValue_noSignal()
	{
		OrbitGeometry geom;
		double origE = geom.eccentricity();
		QSignalSpy spy(&geom, &OrbitGeometry::geometryChanged);
		geom.setEccentricity(origE);
		QCOMPARE(spy.count(), 0);
	}

	// --- Perihelion position ---

	void perihelion_circularOrbit()
	{
		// Circular orbit: e=0, a=1 AU, all angles=0
		OrbitGeometry geom;
		geom.setEccentricity(0.0);
		geom.setSemiMajorAxis(1.0);
		geom.setInclination(0.0);
		geom.setAscendingNode(0.0);
		geom.setArgPericentre(0.0);

		QVector3D peri = geom.perihelionPos();
		float dist = peri.length();
		// Distance from Sun should be 1.0 AU
		QVERIFY2(
		    qAbs(dist - 1.0f) < 0.02f,
		    qPrintable(
		        QString("perihelion dist=%1, pos=(%2,%3,%4)").arg(dist).arg(peri.x()).arg(peri.y()).arg(peri.z())));
	}

	void perihelion_ellipticalOrbit()
	{
		// e=0.5, a=2.0 AU -> q = a*(1-e) = 1.0
		OrbitGeometry geom;
		geom.setEccentricity(0.5);
		geom.setSemiMajorAxis(2.0);
		geom.setInclination(0.0);
		geom.setAscendingNode(0.0);
		geom.setArgPericentre(0.0);

		QVector3D peri = geom.perihelionPos();
		float dist = peri.length();
		QVERIFY(qAbs(dist - 1.0f) < 0.02f);  // q = 1.0 AU
	}

	void perihelion_hyperbolicOrbit()
	{
		// e=3.3769, a=-0.8449 -> q = |a|*(e-1) = 0.8449*2.3769 = 2.008
		OrbitGeometry geom;
		geom.setEccentricity(3.3769);
		geom.setSemiMajorAxis(-0.8449);
		geom.setInclination(44.027);
		geom.setAscendingNode(308.183);
		geom.setArgPericentre(209.103);

		QVector3D peri = geom.perihelionPos();
		float dist = peri.length();
		QVERIFY(qAbs(dist - 2.008f) < 0.05f);  // q ≈ 2.008 AU
	}
};

#endif  // TST_ORBITGEOMETRY_H
