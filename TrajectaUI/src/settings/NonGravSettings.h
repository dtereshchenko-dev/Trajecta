#pragma once

#include <QObject>
#include <QQmlEngine>
#include <QSettings>

// Non-gravitational force model and fitted coefficient selection. Ephemerides always
// use the accepted solution's force model; a gravity-only forecast needs a
// gravity-only fit, not a fitted state propagated without its parameters.
class NonGravSettings : public QObject {
	Q_OBJECT
	QML_ELEMENT

	Q_PROPERTY(bool enabled MEMBER m_enabled NOTIFY changed)
	Q_PROPERTY(int model MEMBER m_model NOTIFY changed)
	Q_PROPERTY(int solveMask MEMBER m_solveMask NOTIFY changed)

   public:
	using QObject::QObject;

	[[nodiscard]] bool enabled() const noexcept { return m_enabled; }

	[[nodiscard]] int model() const noexcept { return m_model; }

	[[nodiscard]] int solveMask() const noexcept { return m_solveMask; }

	void save(QSettings& s) const;
	void load(QSettings& s);

   signals:
	void changed();

   private:
	bool m_enabled{false};
	int m_model{1};      // 1 = MarsdenH2O default
	int m_solveMask{0};  // 0 = model default
};
