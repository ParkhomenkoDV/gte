package logger

import (
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"time"
)

// Level = уровень логирования
type Level int

const (
	NOLOG Level = iota - 1
	STDOUT
	DEBUG
	INFO
	WARN
	ERROR
	FATAL
)

var levelNames = map[Level]string{
	NOLOG:  "NOLOG",
	STDOUT: "STDOUT",
	DEBUG:  "DEBUG",
	INFO:   "INFO",
	WARN:   "WARN",
	ERROR:  "ERROR",
	FATAL:  "FATAL",
}

// Log = запись лога
type Log struct {
	Timestamp time.Time              `json:"timestamp"`
	Level     string                 `json:"level"`
	Message   string                 `json:"message"`
	Fields    map[string]interface{} `json:"fields,omitempty"`
	File      string                 `json:"file,omitempty"`
	Line      int                    `json:"line,omitempty"`
	Function  string                 `json:"function,omitempty"`
}

// Logger основной логгер
type Logger struct {
	level  Level
	output io.Writer
	mu     sync.Mutex
	fields map[string]interface{}
}

// New создает новый логгер
func New(levelName string, outputFile string) (*Logger, error) {
	// Определяем уровень логирования
	var level Level
	switch strings.ToLower(levelName) {
	case "fatal":
		level = FATAL
	case "error":
		level = ERROR
	case "warn", "warning":
		level = WARN
	case "info":
		level = INFO
	case "debug":
		level = DEBUG
	case "stdout":
		return &Logger{level: STDOUT, output: os.Stdout}, nil
	default: // nolog
		return &Logger{level: NOLOG, output: io.Discard}, nil
	}

	// Настраиваем вывод
	output, err := os.OpenFile(outputFile, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0644)
	if err != nil {
		return &Logger{
			level:  level,
			output: os.Stderr,
			fields: make(map[string]interface{}),
		}, fmt.Errorf("открытие файла логов: %v", err)
	}

	return &Logger{
		level:  level,
		output: output,
		fields: make(map[string]interface{}),
	}, nil
}

// SetLevel устанавливает уровень логирования
func (l *Logger) SetLevel(level Level) {
	l.mu.Lock()
	defer l.mu.Unlock()
	l.level = level
}

// SetOutput устанавливает вывод
func (l *Logger) SetOutput(w io.Writer) {
	l.mu.Lock()
	defer l.mu.Unlock()
	l.output = w
}

// WithFields добавляет постоянные поля к логгеру
func (l *Logger) WithFields(fields map[string]interface{}) *Logger {
	l.mu.Lock()
	defer l.mu.Unlock()

	newFields := make(map[string]interface{})
	for k, v := range l.fields {
		newFields[k] = v
	}
	for k, v := range fields {
		newFields[k] = v
	}

	return &Logger{
		level:  l.level,
		output: l.output,
		fields: newFields,
	}
}

// log записывает сообщение
func (l *Logger) log(level Level, msg string, fields map[string]interface{}) {
	if level < l.level {
		return
	}

	// Проверяем, есть ли output
	l.mu.Lock()
	output := l.output
	l.mu.Unlock()
	if output == nil {
		return
	}

	// Получаем информацию о caller
	pc, file, line, ok := runtime.Caller(2)
	var funcName string
	if ok {
		file = filepath.Base(file)
		if fn := runtime.FuncForPC(pc); fn != nil {
			funcName = fn.Name()
			// Оставляем только имя функции
			if idx := strings.LastIndex(funcName, "."); idx != -1 {
				funcName = funcName[idx+1:]
			}
		}
	}

	// Создаем запись
	entry := Log{
		Timestamp: time.Now().UTC(),
		Level:     levelNames[level],
		Message:   msg,
		File:      file,
		Line:      line,
		Function:  funcName,
	}

	// Объединяем поля
	allFields := make(map[string]interface{})
	for k, v := range l.fields {
		allFields[k] = v
	}
	for k, v := range fields {
		allFields[k] = v
	}
	if len(allFields) > 0 {
		entry.Fields = allFields
	}

	l.mu.Lock()
	defer l.mu.Unlock()

	data, err := json.Marshal(entry)
	if err != nil {
		// Если не можем замаршалить в JSON, пишем просто текст
		fmt.Fprintf(l.output, "[%s] %s: %s\n",
			entry.Timestamp.Format(time.RFC3339),
			entry.Level,
			msg)
	} else {
		fmt.Fprintln(l.output, string(data))
	}
}

// Debug логирует отладочное сообщение
func (l *Logger) Debug(msg string, fields ...map[string]interface{}) {
	l.log(DEBUG, msg, mergeFields(fields))
}

// Info логирует информационное сообщение
func (l *Logger) Info(msg string, fields ...map[string]interface{}) {
	l.log(INFO, msg, mergeFields(fields))
}

// Warn логирует предупреждение
func (l *Logger) Warn(msg string, fields ...map[string]interface{}) {
	l.log(WARN, msg, mergeFields(fields))
}

// Error логирует ошибку
func (l *Logger) Error(msg string, fields ...map[string]interface{}) {
	l.log(ERROR, msg, mergeFields(fields))
}

// Fatal логирует фатальную ошибку и завершает программу
func (l *Logger) Fatal(msg string, fields ...map[string]interface{}) {
	l.log(FATAL, msg, mergeFields(fields))
	os.Exit(1)
}

// mergeFields объединяет несколько мап полей
func mergeFields(fields []map[string]interface{}) map[string]interface{} {
	if len(fields) == 0 {
		return nil
	}
	if len(fields) == 1 {
		return fields[0]
	}

	result := make(map[string]interface{})
	for _, f := range fields {
		for k, v := range f {
			result[k] = v
		}
	}
	return result
}
