#pragma once

#include <QObject>
#include <QSettings>
#include <QQmlEngine>

class SigmaClipSettings : public QObject {
	Q_OBJECT
	QML_ELEMENT

	Q_PROPERTY(bool enabled MEMBER m_enabled NOTIFY changed)
	Q_PROPERTY(double threshold MEMBER m_threshold NOTIFY changed)

   public:
	using QObject::QObject;

	[[nodiscard]] bool enabled() const noexcept { return m_enabled; }

	[[nodiscard]] double threshold() const noexcept { return m_threshold; }

	void save(QSettings& s) const;
	void load(QSettings& s);

   signals:
	void changed();

   private:
	bool m_enabled{true};
	double m_threshold{3.0};
};
