#ifndef TST_GRIDGEOMETRY_H
#define TST_GRIDGEOMETRY_H

#include <QTest>
#include <QSignalSpy>
#include "geometry/GridGeometry.h"

class TestGridGeometry : public QObject {
	Q_OBJECT

   private slots:

	void defaultValues()
	{
		GridGeometry geom;
		QCOMPARE(geom.size(), 10.0f);
		QCOMPARE(geom.divisions(), 20);
	}

	void setSize_emitsSignal()
	{
		GridGeometry geom;
		QSignalSpy spy(&geom, &GridGeometry::geometryChanged);
		geom.setSize(5.0f);
		QCOMPARE(geom.size(), 5.0f);
		QVERIFY(spy.count() >= 1);
	}

	void setDivisions_emitsSignal()
	{
		GridGeometry geom;
		QSignalSpy spy(&geom, &GridGeometry::geometryChanged);
		geom.setDivisions(10);
		QCOMPARE(geom.divisions(), 10);
		QVERIFY(spy.count() >= 1);
	}

	void setSameSize_noSignal()
	{
		GridGeometry geom;
		QSignalSpy spy(&geom, &GridGeometry::geometryChanged);
		geom.setSize(10.0f);  // default
		QCOMPARE(spy.count(), 0);
	}

	void setSameDivisions_noSignal()
	{
		GridGeometry geom;
		QSignalSpy spy(&geom, &GridGeometry::geometryChanged);
		geom.setDivisions(20);  // default
		QCOMPARE(spy.count(), 0);
	}
};

#endif  // TST_GRIDGEOMETRY_H
