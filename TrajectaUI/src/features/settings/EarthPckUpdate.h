#pragma once

#include <QObject>
#include <QQmlEngine>
#include <QNetworkAccessManager>
#include <QPointer>
#include <QUrl>
#include <QVariantMap>
#include <array>

class QNetworkReply;

class EarthPckUpdate : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Session-owned")
	Q_PROPERTY(QVariantMap states READ states NOTIFY changed)
	Q_PROPERTY(int availableCount READ availableCount NOTIFY changed)

   public:
	explicit EarthPckUpdate(QObject* parent = nullptr);
	EarthPckUpdate(QUrl directoryUrl, QString directory, QString metadataPath, QObject* parent = nullptr);
	static void applyPending(const QString& directory, const QString& metadataPath);
	static void applyPending();
	QVariantMap states() const;
	int availableCount() const;
	Q_INVOKABLE void check();
	Q_INVOKABLE void install(const QString& id);

   signals:
	void changed();
	void installed();

   private:
	struct Source {
		QString id;
		QString file;
		QString remote;
		QString originalHash;
		QString message;
		bool available = false;
		bool checking = false;
		bool downloading = false;
	};

	Source* source(const QString& id);
	QNetworkReply* get(const QUrl& url);
	void finishCheck(QNetworkReply* reply);
	void finishDownload(QNetworkReply* reply, Source* source);
	QNetworkAccessManager m_network;
	QUrl m_directoryUrl;
	QString m_directory;
	QString m_metadataPath;
	QPointer<QNetworkReply> m_check;
	std::array<Source, 3> m_sources;
};
