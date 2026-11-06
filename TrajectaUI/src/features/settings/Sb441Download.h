#pragma once

#include <QCryptographicHash>
#include <QNetworkAccessManager>
#include <QObject>
#include <QQmlEngine>
#include <QUrl>
#include <memory>

class QNetworkReply;
class QSaveFile;

class Sb441Download : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Session-owned")
	Q_PROPERTY(bool running READ running NOTIFY changed)
	Q_PROPERTY(double progress READ progress NOTIFY changed)
	Q_PROPERTY(QString message READ message NOTIFY changed)

   public:
	explicit Sb441Download(QObject* parent = nullptr);
	Sb441Download(QUrl source, QString directory, qint64 expectedSize, QByteArray expectedSha256,
	              QObject* parent = nullptr);
	~Sb441Download() override;
	[[nodiscard]] static QString installedPath();

	[[nodiscard]] bool running() const { return m_reply != nullptr; }

	[[nodiscard]] double progress() const;

	[[nodiscard]] QString message() const { return m_message; }

	Q_INVOKABLE void start();
	Q_INVOKABLE void cancel();

   signals:
	void changed();
	void installed();

   private:
	void consume();
	void finish();
	QNetworkAccessManager m_network;
	QNetworkReply* m_reply = nullptr;
	std::unique_ptr<QSaveFile> m_file;
	QCryptographicHash m_hash{QCryptographicHash::Sha256};
	QUrl m_source;
	QString m_directory;
	qint64 m_expectedSize;
	QByteArray m_expectedSha256;
	qint64 m_received = 0;
	QString m_message;
	QString m_failure;
	bool m_cancelled = false;
	bool m_finishing = false;
};
