#pragma once

#include <QObject>
#include <QProcess>
#include <QQmlEngine>
#include <QTimer>

class ApplicationUpdate : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Session-owned")
	Q_PROPERTY(bool available READ available NOTIFY changed)
	Q_PROPERTY(QString version READ version NOTIFY changed)
	Q_PROPERTY(QString error READ error NOTIFY changed)

   public:
	explicit ApplicationUpdate(QObject* parent = nullptr);
	~ApplicationUpdate() override;
	static QString availableVersion(const QByteArray& output, const QString& currentVersion);

	bool available() const { return !m_version.isEmpty(); }

	QString version() const { return m_version; }

	QString error() const { return m_error; }

	Q_INVOKABLE void check();
	Q_INVOKABLE bool openUpdater();

   signals:
	void changed();

   private:
	QString toolPath() const;
	QProcess m_process;
	QTimer m_timeout;
	QString m_version;
	QString m_error;
};
