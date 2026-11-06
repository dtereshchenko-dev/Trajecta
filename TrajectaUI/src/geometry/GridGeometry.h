#ifndef GRIDGEOMETRY_H
#define GRIDGEOMETRY_H

#include <QQuick3DGeometry>
#include <QQmlEngine>

class GridGeometry : public QQuick3DGeometry {
	Q_OBJECT
	QML_ELEMENT

	Q_PROPERTY(float size READ size WRITE setSize NOTIFY geometryChanged)
	Q_PROPERTY(int divisions READ divisions WRITE setDivisions NOTIFY geometryChanged)

   public:
	explicit GridGeometry(QQuick3DObject* parent = nullptr);

	float size() const { return m_size; }

	int divisions() const { return m_divisions; }

	void setSize(float v)
	{
		if (m_size != v) {
			m_size = v;
			rebuild();
			emit geometryChanged();
		}
	}

	void setDivisions(int v)
	{
		if (m_divisions != v) {
			m_divisions = v;
			rebuild();
			emit geometryChanged();
		}
	}

   private:
	void rebuild();
	float m_size = 10.0f;
	int m_divisions = 20;
};

#endif  // GRIDGEOMETRY_H
