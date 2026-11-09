#pragma once

#include <QFutureWatcher>
#include <QObject>
#include <QQmlEngine>
#include <QString>
#include <QVariantMap>
#include "orbital_core/obs/ObservatoryUpdates.h"
#include <functional>
#include <stop_token>

class ObservatoryUpdate : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Session-owned")
	Q_PROPERTY(int availableCount READ availableCount NOTIFY changed)
	Q_PROPERTY(QVariantMap states READ states NOTIFY changed)

   public:
	explicit ObservatoryUpdate(QObject* parent = nullptr);
	using Fetch =
	    std::function<orbital_core::obs::ObservatoryUpdate(const QString&, std::string_view, std::stop_token)>;
	explicit ObservatoryUpdate(Fetch fetch, QObject* parent = nullptr);
	~ObservatoryUpdate() override;
	int availableCount() const;
	QVariantMap states() const;
	Q_INVOKABLE void check();
	Q_INVOKABLE void install(const QString& id);

   signals:
	void changed();
	void installed();

   private:
	struct Result {
		QByteArray original;
		QByteArray contents;
		QString error;
		std::size_t updated{};
		std::size_t unmapped{};
		bool existed = false;
		bool changed = false;
	};

	struct Source {
		QString id;
		QFutureWatcher<Result> watcher;
		std::stop_source cancel;
		Result result;
		QString message;
		bool available = false;
		bool checking = false;
	};

	Source* source(const QString& id);
	void checkSource(Source& source);
	Fetch m_fetch;
	Source m_optical;
	Source m_radar;
};
