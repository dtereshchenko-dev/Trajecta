#pragma once

#include <QQuick3DGeometry>
#include <QQmlEngine>
#include <QVariantList>

class PolylineGeometry : public QQuick3DGeometry {
	Q_OBJECT
	QML_ELEMENT
	Q_PROPERTY(QVariantList points READ points WRITE setPoints NOTIFY pointsChanged)
	Q_PROPERTY(bool independentSegments READ independentSegments WRITE setIndependentSegments NOTIFY
	               independentSegmentsChanged)

   public:
	explicit PolylineGeometry(QQuick3DObject* parent = nullptr);

	QVariantList points() const { return m_points; }

	void setPoints(const QVariantList& points);

	bool independentSegments() const { return m_independentSegments; }

	void setIndependentSegments(bool value);

   signals:
	void pointsChanged();
	void independentSegmentsChanged();

   private:
	QVariantList m_points;
	bool m_independentSegments{};
	void rebuild();
};
