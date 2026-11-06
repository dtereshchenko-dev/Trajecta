#ifndef ORBITGEOMETRY_H
#define ORBITGEOMETRY_H

#include <QQuick3DGeometry>
#include <QQmlEngine>
#include <QVector3D>

class OrbitGeometry : public QQuick3DGeometry {
	Q_OBJECT
	QML_ELEMENT

	Q_PROPERTY(double eccentricity READ eccentricity WRITE setEccentricity NOTIFY geometryChanged)
	Q_PROPERTY(double semiMajorAxis READ semiMajorAxis WRITE setSemiMajorAxis NOTIFY geometryChanged)
	Q_PROPERTY(double perigeeAu READ perigeeAu WRITE setPerigeeAu NOTIFY geometryChanged)
	Q_PROPERTY(double inclination READ inclination WRITE setInclination NOTIFY geometryChanged)
	Q_PROPERTY(double ascendingNode READ ascendingNode WRITE setAscendingNode NOTIFY geometryChanged)
	Q_PROPERTY(double argPericentre READ argPericentre WRITE setArgPericentre NOTIFY geometryChanged)
	Q_PROPERTY(int segments READ segments WRITE setSegments NOTIFY geometryChanged)
	Q_PROPERTY(QVector3D perihelionPos READ perihelionPos NOTIFY geometryChanged)

   public:
	explicit OrbitGeometry(QQuick3DObject* parent = nullptr);

	double eccentricity() const { return m_e; }

	double semiMajorAxis() const { return m_a; }

	double perigeeAu() const { return m_q; }

	double inclination() const { return m_i; }

	double ascendingNode() const { return m_Omega; }

	double argPericentre() const { return m_omega; }

	int segments() const { return m_segments; }

	void setEccentricity(double v)
	{
		if (m_e != v) {
			m_e = v;
			rebuild();
			emit geometryChanged();
		}
	}

	void setSemiMajorAxis(double v)
	{
		if (m_a != v) {
			m_a = v;
			rebuild();
			emit geometryChanged();
		}
	}

	void setPerigeeAu(double v)
	{
		if (m_q != v) {
			m_q = v;
			rebuild();
			emit geometryChanged();
		}
	}

	void setInclination(double v)
	{
		if (m_i != v) {
			m_i = v;
			rebuild();
			emit geometryChanged();
		}
	}

	void setAscendingNode(double v)
	{
		if (m_Omega != v) {
			m_Omega = v;
			rebuild();
			emit geometryChanged();
		}
	}

	void setArgPericentre(double v)
	{
		if (m_omega != v) {
			m_omega = v;
			rebuild();
			emit geometryChanged();
		}
	}

	void setSegments(int v)
	{
		if (m_segments != v) {
			m_segments = v;
			rebuild();
			emit geometryChanged();
		}
	}

	QVector3D perihelionPos() const { return m_perihelionPos; }

   private:
	void rebuild();
	QVector3D orbitalToEcliptic(double x_orb, double y_orb) const;
	QVector3D m_perihelionPos;

	double m_e = 3.3769;
	double m_a = -0.8449;
	double m_q = 0.0;          // AU; zero uses the conventional (a, e) geometry
	double m_i = 44.027;       // degrees
	double m_Omega = 308.183;  // degrees
	double m_omega = 209.103;  // degrees
	int m_segments = 200;
};

#endif  // ORBITGEOMETRY_H
