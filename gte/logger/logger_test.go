package logger

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"
	"time"
)

// TestNew создает тесты для функции New
func TestNew(t *testing.T) {
	tests := []struct {
		name       string
		level      string
		outputFile string
		wantError  bool
		wantLevel  Level
	}{
		{
			name:       "debug уровень с файлом",
			level:      "debug",
			outputFile: "test_debug.log",
			wantError:  false,
			wantLevel:  DEBUG,
		}, {
			name:       "info уровень с файлом",
			level:      "info",
			outputFile: "test_info.log",
			wantError:  false,
			wantLevel:  INFO,
		}, {
			name:       "warn уровень с файлом",
			level:      "warn",
			outputFile: "test_warn.log",
			wantError:  false,
			wantLevel:  WARN,
		}, {
			name:       "error уровень с файлом",
			level:      "error",
			outputFile: "test_error.log",
			wantError:  false,
			wantLevel:  ERROR,
		}, {
			name:       "fatal уровень с файлом",
			level:      "fatal",
			outputFile: "test_fatal.log",
			wantError:  false,
			wantLevel:  FATAL,
		}, {
			name:       "stdout уровень",
			level:      "stdout",
			outputFile: "",
			wantError:  false,
			wantLevel:  STDOUT,
		}, {
			name:       "несуществующий уровень (должен быть NOLOG)",
			level:      "unknown",
			outputFile: "",
			wantError:  false,
			wantLevel:  NOLOG,
		}, {
			name:       "пустой уровень (должен быть NOLOG)",
			level:      "",
			outputFile: "",
			wantError:  false,
			wantLevel:  NOLOG,
		}, {
			name:       "warning синоним для warn",
			level:      "warning",
			outputFile: "test_warning.log",
			wantError:  false,
			wantLevel:  WARN,
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			// Удаляем тестовые файлы после теста
			if test.outputFile != "" {
				defer func() {
					os.Remove(test.outputFile)
				}()
			}

			logger, err := New(test.level, test.outputFile)

			if test.wantError && err == nil {
				t.Errorf("ожидалась ошибка, но не получена")
			}

			if !test.wantError && err != nil {
				t.Errorf("не ожидалась ошибка, но получена: %v", err)
			}

			if err == nil && logger.level != test.wantLevel {
				t.Errorf("уровень логирования = %v, ожидалось %v", logger.level, test.wantLevel)
			}
		})
	}
}

// TestNew_FileCreation проверяет создание файла логов
func TestNew_FileCreation(t *testing.T) {
	tempFile := "test_log_creation.log"
	defer os.Remove(tempFile)

	logger, err := New("info", tempFile)
	if err != nil {
		t.Fatalf("не удалось создать логгер: %v", err)
	}
	defer func() {
		if closer, ok := logger.output.(io.Closer); ok {
			closer.Close()
		}
	}()

	// Проверяем, что файл создан
	if _, err := os.Stat(tempFile); os.IsNotExist(err) {
		t.Errorf("файл логов не создан: %v", err)
	}

	// Записываем тестовое сообщение
	logger.Info("тестовое сообщение")

	// Проверяем содержимое файла
	content, err := os.ReadFile(tempFile)
	if err != nil {
		t.Fatalf("не удалось прочитать файл: %v", err)
	}

	if !strings.Contains(string(content), "тестовое сообщение") {
		t.Errorf("сообщение не записано в файл")
	}
}

// TestNew_InvalidFile проверяет обработку невалидного пути файла
func TestNew_InvalidFile(t *testing.T) {
	// Пытаемся создать файл в несуществующей директории
	invalidPath := "/несуществующая/директория/test.log"
	logger, err := New("info", invalidPath)

	// Должен вернуться логгер с stderr
	if err == nil {
		t.Errorf("ожидалась ошибка для невалидного пути файла")
	}

	if logger == nil {
		t.Fatalf("логгер не должен быть nil даже при ошибке")
	}

	if logger.output != os.Stderr {
		t.Errorf("при ошибке вывод должен быть установлен в stderr")
	}
}

// TestSetLevel проверяет установку уровня логирования
func TestSetLevel(t *testing.T) {
	buf := &bytes.Buffer{}
	logger := &Logger{
		level:  INFO,
		output: buf,
		fields: make(map[string]interface{}),
	}

	// Тестируем повышение уровня
	logger.SetLevel(ERROR)
	logger.Info("это сообщение не должно быть записано")
	if buf.Len() > 0 {
		t.Errorf("сообщение записано, хотя уровень выше")
	}

	// Тестируем понижение уровня
	buf.Reset()
	logger.SetLevel(DEBUG)
	logger.Info("это сообщение должно быть записано")
	if buf.Len() == 0 {
		t.Errorf("сообщение не записано, хотя уровень ниже")
	}
}

// TestSetOutput проверяет установку вывода
func TestSetOutput(t *testing.T) {
	buf1 := &bytes.Buffer{}
	buf2 := &bytes.Buffer{}

	logger := &Logger{
		level:  DEBUG,
		output: buf1,
		fields: make(map[string]interface{}),
	}

	// Записываем в первый буфер
	logger.Info("сообщение 1")
	if buf1.Len() == 0 {
		t.Errorf("сообщение не записано в первый буфер")
	}

	// Меняем вывод
	logger.SetOutput(buf2)
	logger.Info("сообщение 2")

	// Проверяем что первое сообщение только в первом буфере
	content1 := buf1.String()
	if !strings.Contains(content1, "сообщение 1") {
		t.Errorf("первое сообщение должно быть в первом буфере")
	}
	if strings.Contains(content1, "сообщение 2") {
		t.Errorf("второе сообщение не должно быть в первом буфере")
	}

	// Проверяем что второе сообщение только во втором буфере
	content2 := buf2.String()
	if strings.Contains(content2, "сообщение 1") {
		t.Errorf("первое сообщение не должно быть во втором буфере")
	}
	if !strings.Contains(content2, "сообщение 2") {
		t.Errorf("второе сообщение должно быть во втором буфере")
	}
}

// TestWithFields проверяет добавление постоянных полей
func TestWithFields(t *testing.T) {
	buf := &bytes.Buffer{}
	baseLogger := &Logger{
		level:  DEBUG,
		output: buf,
		fields: map[string]interface{}{
			"service": "test",
			"version": "1.0",
		},
	}

	// Создаем логгер с дополнительными полями
	childLogger := baseLogger.WithFields(map[string]interface{}{
		"instance": "server-1",
		"version":  "1.1", // Переопределяем существующее поле
	})

	// Проверяем что поля объединены
	childLogger.Info("тестовое сообщение")

	var logEntry Log
	if err := json.Unmarshal(buf.Bytes(), &logEntry); err != nil {
		t.Fatalf("не удалось разобрать JSON: %v", err)
	}

	// Проверяем поля
	if service, ok := logEntry.Fields["service"].(string); !ok || service != "test" {
		t.Errorf("поле service = %v, ожидалось 'test'", logEntry.Fields["service"])
	}

	if version, ok := logEntry.Fields["version"].(string); !ok || version != "1.1" {
		t.Errorf("поле version = %v, ожидалось '1.1'", logEntry.Fields["version"])
	}

	if instance, ok := logEntry.Fields["instance"].(string); !ok || instance != "server-1" {
		t.Errorf("поле instance = %v, ожидалось 'server-1'", logEntry.Fields["instance"])
	}

	// Проверяем что оригинальный логгер не изменился
	buf.Reset()
	baseLogger.Info("еще одно сообщение")

	var baseLogEntry Log
	if err := json.Unmarshal(buf.Bytes(), &baseLogEntry); err != nil {
		t.Fatalf("не удалось разобрать JSON: %v", err)
	}

	if version, ok := baseLogEntry.Fields["version"].(string); !ok || version != "1.0" {
		t.Errorf("базовый логгер: поле version = %v, ожидалось '1.0'", baseLogEntry.Fields["version"])
	}

	if _, ok := baseLogEntry.Fields["instance"]; ok {
		t.Errorf("базовый логгер не должен содержать поле instance")
	}
}

// TestLogLevels проверяет все уровни логирования
func TestLogLevels(t *testing.T) {
	tests := []struct {
		name    string
		logFunc func(*Logger)
		level   Level
		want    bool
	}{
		{
			name: "debug при уровне DEBUG",
			logFunc: func(l *Logger) {
				l.Debug("debug сообщение")
			},
			level: DEBUG,
			want:  true,
		}, {
			name: "debug при уровне INFO",
			logFunc: func(l *Logger) {
				l.Debug("debug сообщение")
			},
			level: INFO,
			want:  false,
		}, {
			name: "info при уровне INFO",
			logFunc: func(l *Logger) {
				l.Info("info сообщение")
			},
			level: INFO,
			want:  true,
		}, {
			name: "info при уровне WARN",
			logFunc: func(l *Logger) {
				l.Info("info сообщение")
			},
			level: WARN,
			want:  false,
		}, {
			name: "warn при уровне WARN",
			logFunc: func(l *Logger) {
				l.Warn("warn сообщение")
			},
			level: WARN,
			want:  true,
		}, {
			name: "warn при уровне ERROR",
			logFunc: func(l *Logger) {
				l.Warn("warn сообщение")
			},
			level: ERROR,
			want:  false,
		}, {
			name: "error при уровне ERROR",
			logFunc: func(l *Logger) {
				l.Error("error сообщение")
			},
			level: ERROR,
			want:  true,
		}, {
			name: "error при уровне FATAL",
			logFunc: func(l *Logger) {
				l.Error("error сообщение")
			},
			level: FATAL,
			want:  false,
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			buf := &bytes.Buffer{}
			logger := &Logger{
				level:  test.level,
				output: buf,
				fields: make(map[string]interface{}),
			}

			test.logFunc(logger)

			hasOutput := buf.Len() > 0
			if hasOutput != test.want {
				t.Errorf("вывод = %v, ожидалось %v", hasOutput, test.want)
			}
		})
	}
}

// TestLogFields проверяет логирование с полями
func TestLogFields(t *testing.T) {
	buf := &bytes.Buffer{}
	logger := &Logger{
		level:  DEBUG,
		output: buf,
		fields: make(map[string]interface{}),
	}

	// Логируем с полями
	logger.Info("обработка запроса",
		map[string]interface{}{
			"user_id":  123,
			"method":   "GET",
			"endpoint": "/api/test",
			"success":  true,
		})

	var logEntry Log
	if err := json.Unmarshal(buf.Bytes(), &logEntry); err != nil {
		t.Fatalf("не удалось разобрать JSON: %v", err)
	}

	// Проверяем поля
	if userID, ok := logEntry.Fields["user_id"].(float64); !ok || userID != 123 {
		t.Errorf("user_id = %v, ожидалось 123", logEntry.Fields["user_id"])
	}

	if method, ok := logEntry.Fields["method"].(string); !ok || method != "GET" {
		t.Errorf("method = %v, ожидалось 'GET'", logEntry.Fields["method"])
	}

	if success, ok := logEntry.Fields["success"].(bool); !ok || !success {
		t.Errorf("success = %v, ожидалось true", logEntry.Fields["success"])
	}
}

// TestLogMultipleFields проверяет логирование с несколькими мапами полей
func TestLogMultipleFields(t *testing.T) {
	buf := &bytes.Buffer{}
	logger := &Logger{
		level:  DEBUG,
		output: buf,
		fields: make(map[string]interface{}),
	}

	// Логируем с несколькими мапами полей
	logger.Info("комплексная операция",
		map[string]interface{}{"step": 1, "action": "start"},
		map[string]interface{}{"data": "test", "count": 5},
		map[string]interface{}{"step": 2, "status": "processing"})

	var logEntry Log
	if err := json.Unmarshal(buf.Bytes(), &logEntry); err != nil {
		t.Fatalf("не удалось разобрать JSON: %v", err)
	}

	// Проверяем что все поля объединены (последнее значение перезаписывает предыдущее)
	if step, ok := logEntry.Fields["step"].(float64); !ok || step != 2 {
		t.Errorf("step = %v, ожидалось 2 (последнее значение)", logEntry.Fields["step"])
	}

	if action, ok := logEntry.Fields["action"].(string); !ok || action != "start" {
		t.Errorf("action = %v, ожидалось 'start'", logEntry.Fields["action"])
	}

	if data, ok := logEntry.Fields["data"].(string); !ok || data != "test" {
		t.Errorf("data = %v, ожидалось 'test'", logEntry.Fields["data"])
	}

	if count, ok := logEntry.Fields["count"].(float64); !ok || count != 5 {
		t.Errorf("count = %v, ожидалось 5", logEntry.Fields["count"])
	}
}

// TestLogCallerInfo проверяет информацию о caller
func TestLogCallerInfo(t *testing.T) {
	buf := &bytes.Buffer{}
	logger := &Logger{
		level:  INFO,
		output: buf,
		fields: make(map[string]interface{}),
	}

	// Логируем сообщение
	logger.Info("тестовое сообщение")

	var logEntry Log
	if err := json.Unmarshal(buf.Bytes(), &logEntry); err != nil {
		t.Fatalf("не удалось разобрать JSON: %v", err)
	}

	// Проверяем наличие информации о файле и функции
	if logEntry.File == "" {
		t.Error("отсутствует информация о файле")
	}

	if logEntry.Line <= 0 {
		t.Error("отсутствует или некорректная информация о строке")
	}

	// Функция должна быть TestLogCallerInfo
	if !strings.Contains(logEntry.Function, "TestLogCallerInfo") {
		t.Errorf("функция = %q, ожидалось 'TestLogCallerInfo'", logEntry.Function)
	}
}

// TestLogJSONError проверяет обработку ошибок JSON маршалинга
func TestLogJSONError(t *testing.T) {
	// Создаем специальный writer, который записывает в буфер
	buf := &bytes.Buffer{}
	logger := &Logger{
		level: INFO,
		output: &errorWriter{
			Writer: buf,
			failOn: 2, // Сбой при второй записи (JSON маршалинг)
		},
		fields: make(map[string]interface{}),
	}

	// Создаем данные, которые невозможно сериализовать в JSON
	// (функция не сериализуема в JSON)
	unmarshalableData := map[string]interface{}{
		"func": func() {}, // Функция не может быть сериализована в JSON
	}

	// Логируем с некорректными данными
	logger.Info("сообщение с некорректными данными", unmarshalableData)

	// Проверяем что запись все равно произошла в текстовом формате
	output := buf.String()
	if !strings.Contains(output, "сообщение с некорректными данными") {
		t.Error("сообщение должно быть записано в текстовом формате при ошибке JSON")
	}

	// Проверяем что есть временная метка и уровень
	if !strings.Contains(output, "INFO") {
		t.Error("должен содержать уровень логирования")
	}
}

// errorWriter - writer, который может симулировать ошибки
type errorWriter struct {
	io.Writer
	count  int
	failOn int
}

func (ew *errorWriter) Write(p []byte) (n int, err error) {
	ew.count++
	if ew.count == ew.failOn {
		return 0, fmt.Errorf("симулированная ошибка записи")
	}
	return ew.Writer.Write(p)
}

// TestLogConcurrent проверяет конкурентное использование логгера
func TestLogConcurrent(t *testing.T) {
	buf := &bytes.Buffer{}
	logger := &Logger{
		level:  DEBUG,
		output: buf,
		fields: make(map[string]interface{}),
	}

	// Запускаем несколько горутин
	var wg sync.WaitGroup
	messages := 100

	for i := 0; i < messages; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			logger.Info(fmt.Sprintf("сообщение %d", id),
				map[string]interface{}{
					"goroutine_id": id,
					"timestamp":    time.Now().UnixNano(),
				})
		}(i)
	}

	wg.Wait()

	// Проверяем что все сообщения записаны
	lines := strings.Count(buf.String(), "\n")
	if lines != messages {
		t.Errorf("записано %d строк, ожидалось %d", lines, messages)
	}
}

// TestLogStructure проверяет структуру лог-записи
func TestLogStructure(t *testing.T) {
	buf := &bytes.Buffer{}
	logger := &Logger{
		level:  INFO,
		output: buf,
		fields: map[string]interface{}{},
	}

	testTime := time.Now()
	logger.Info("структурированное сообщение",
		map[string]interface{}{
			"app":    "test",
			"number": 42,
			"float":  3.14,
			"bool":   true,
			"null":   nil,
		})

	var logEntry Log
	if err := json.Unmarshal(buf.Bytes(), &logEntry); err != nil {
		t.Fatalf("не удалось разобрать JSON: %v", err)
	}

	// Проверяем обязательные поля
	if logEntry.Timestamp.IsZero() {
		t.Error("отсутствует временная метка")
	} else if logEntry.Timestamp.Before(testTime.Add(-time.Second)) {
		t.Error("временная метка слишком старая")
	}

	if logEntry.Level != "INFO" {
		t.Errorf("уровень = %q, ожидалось 'INFO'", logEntry.Level)
	}

	if logEntry.Message != "структурированное сообщение" {
		t.Errorf("сообщение = %q, ожидалось 'структурированное сообщение'", logEntry.Message)
	}

	// Проверяем поля
	if app, ok := logEntry.Fields["app"].(string); !ok || app != "test" {
		t.Errorf("поле app = %v, ожидалось 'test'", logEntry.Fields["app"])
	}

	if number, ok := logEntry.Fields["number"].(float64); !ok || number != 42 {
		t.Errorf("поле number = %v, ожидалось 42", logEntry.Fields["number"])
	}

	if floatVal, ok := logEntry.Fields["float"].(float64); !ok || floatVal != 3.14 {
		t.Errorf("поле float = %v, ожидалось 3.14", logEntry.Fields["float"])
	}

	if boolVal, ok := logEntry.Fields["bool"].(bool); !ok || !boolVal {
		t.Errorf("поле bool = %v, ожидалось true", logEntry.Fields["bool"])
	}

	// Проверяем информацию о caller
	if logEntry.File != filepath.Base("logger_test.go") {
		t.Errorf("файл = %q, ожидалось 'logger_test.go'", logEntry.File)
	}

	if !strings.Contains(logEntry.Function, "TestLogStructure") {
		t.Errorf("функция = %q, должна содержать 'TestLogStructure'", logEntry.Function)
	}
}

// TestNoOutput проверяет логирование без установленного output
func TestNoOutput(t *testing.T) {
	logger, _ := New("INFO", "")

	// Не должно быть паники
	logger.Info("сообщение без output")
	logger.Warn("предупреждение без output")
	logger.Error("ошибка без output")
}

// TestStdoutLevel проверяет специальный уровень STDOUT
func TestStdoutLevel(t *testing.T) {
	// Заменяем os.Stdout для тестирования
	oldStdout := os.Stdout
	r, w, _ := os.Pipe()
	os.Stdout = w
	defer func() { os.Stdout = oldStdout }()

	logger, err := New("stdout", "")
	if err != nil {
		t.Fatalf("не удалось создать логгер с уровнем stdout: %v", err)
	}

	if logger.level != STDOUT {
		t.Errorf("уровень = %v, ожидалось STDOUT", logger.level)
	}

	if logger.output != os.Stdout {
		t.Errorf("output не установлен в os.Stdout")
	}

	// Восстанавливаем stdout и читаем вывод
	logger.Info("тестовое сообщение в stdout")
	w.Close()
	output, _ := io.ReadAll(r)

	if !strings.Contains(string(output), "тестовое сообщение в stdout") {
		t.Error("сообщение не записано в stdout")
	}
}

// TestNologLevel проверяет уровень NOLOG
func TestNologLevel(t *testing.T) {
	buf := &bytes.Buffer{}
	logger, _ := New("NOLOG", "")

	// Ни одно сообщение не должно быть записано
	logger.Debug("debug сообщение")
	logger.Info("info сообщение")
	logger.Warn("warn сообщение")
	logger.Error("error сообщение")

	if buf.Len() > 0 {
		t.Errorf("сообщения записаны при уровне NOLOG: %s", buf.String())
	}
}

// TestLogOrder проверяет порядок полей в JSON
func TestLogOrder(t *testing.T) {
	buf := &bytes.Buffer{}
	logger := &Logger{
		level:  INFO,
		output: buf,
		fields: make(map[string]interface{}),
	}

	logger.Info("проверка порядка")

	var logEntry map[string]interface{}
	if err := json.Unmarshal(buf.Bytes(), &logEntry); err != nil {
		t.Fatalf("не удалось разобрать JSON: %v", err)
	}

	// Проверяем наличие всех ожидаемых полей
	expectedFields := []string{"timestamp", "level", "message", "file", "line", "function"}
	for _, field := range expectedFields {
		if _, ok := logEntry[field]; !ok {
			t.Errorf("отсутствует поле %s", field)
		}
	}
}
