#ifndef RECENTFILES_H
#define RECENTFILES_H

#include <QObject>
#include <QStringList>
#include <QUrl>
#include <QQmlEngine>

class RecentFiles : public QObject {
	Q_OBJECT
	QML_ELEMENT
	Q_PROPERTY(QStringList files READ files NOTIFY filesChanged)
	Q_PROPERTY(int count READ count NOTIFY filesChanged)

   public:
	explicit RecentFiles(QObject* parent = nullptr);
	~RecentFiles() = default;

	QStringList files() const;
	int count() const;

	Q_INVOKABLE void addFile(const QUrl& url);
	Q_INVOKABLE void removeFile(const QString& path);
	Q_INVOKABLE void clear();

	Q_INVOKABLE QString fileName(int index) const;
	Q_INVOKABLE QString filePath(int index) const;
	Q_INVOKABLE QString displayPath(int index) const;

   signals:
	void filesChanged();

   private:
	void load();
	void save();

	QStringList m_files;
	static constexpr int MaxRecent = 6;
};

#endif  // RECENTFILES_H
