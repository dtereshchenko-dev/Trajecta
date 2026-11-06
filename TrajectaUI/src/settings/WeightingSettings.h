#pragma once

#include <QObject>
#include <QQmlEngine>
#include <QSettings>

class WeightingSettings : public QObject {
	Q_OBJECT
	QML_ELEMENT

	Q_PROPERTY(bool useWeights MEMBER m_useWeights NOTIFY changed)
	Q_PROPERTY(bool observerSigmas MEMBER m_observerSigmas NOTIFY changed)
	Q_PROPERTY(bool timingEnabled MEMBER m_timingEnabled NOTIFY changed)
	Q_PROPERTY(bool debiasEnabled MEMBER m_debiasEnabled NOTIFY changed)
	Q_PROPERTY(int weightProfile MEMBER m_weightProfile NOTIFY changed)

   public:
	using QObject::QObject;

	[[nodiscard]] bool useWeights() const noexcept { return m_useWeights; }

	[[nodiscard]] bool observerSigmas() const noexcept { return m_observerSigmas; }

	[[nodiscard]] bool timingEnabled() const noexcept { return m_timingEnabled; }

	[[nodiscard]] bool debiasEnabled() const noexcept { return m_debiasEnabled; }

	[[nodiscard]] int weightProfile() const noexcept { return m_weightProfile; }

	void save(QSettings& s) const;
	void load(QSettings& s);

   signals:
	void changed();

   private:
	bool m_useWeights{true};
	bool m_observerSigmas{true};
	// Whether the timing term applies at all.  Deliberately a flag and not a
	// number: the value used for an observation that states no timing
	// uncertainty lives in kDefaultTimingSigmaSeconds, so revising that
	// constant takes effect immediately instead of being shadowed by a stale
	// number in QSettings.  Sweeping the value is a CLI concern.
	bool m_timingEnabled{false};
	bool m_debiasEnabled{true};
	int m_weightProfile{1};
};
