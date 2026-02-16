# Makefile for Python project

# Configuration
PROJECT_NAME = gte
PYTHON = python3
PIP = pip3
VENV_DIR = .venv
VENV_ACTIVATE = $(VENV_DIR)/bin/activate
PYTHON_PATH = $(VENV_DIR)/bin/python
PIP_PATH = $(VENV_DIR)/bin/pip
TEST_DIR = gte
BENCH_DIR = gte
SRC_DIR = gte
REQUIREMENTS = requirements.txt

# Colors
RED    = \033[0;31m
GREEN  = \033[0;32m
YELLOW = \033[0;33m
BLUE   = \033[0;34m
RESET  = \033[0m

# Targets
.PHONY: help venv activate install test bench lint format clean

help:
	@echo "Available commands:"
	@echo "  make venv           - Create virtual environment"
	@echo "  make activate       - Activate virtual environment (prints command)"
	@echo "  make install        - Install production dependeRESETies"
	@echo "  make test           - Run tests"
	@echo "  make lint           - Run linters (flake8, pylint)"
	@echo "  make format         - Format code (black, isort)"
	@echo "  make clean          - Clean project"

venv:
	@echo "$(BLUE)Creating virtual environment...$(RESET)"
	$(PYTHON) -m venv $(VENV_DIR)
	@echo "$(BLUE)Virtual environment created in $(VENV_DIR)$(RESET)"
	@echo "To activate, run:"
	@echo "  source $(VENV_ACTIVATE)"

activate:
	@echo "$(BLUE)Run this command to activate virtual environment:$(RESET)"
	@echo "source $(VENV_ACTIVATE)"

install:
	@echo "$(BLUE)Installing production dependencies...$(RESET)"
	$(PIP_PATH) install --upgrade -r $(REQUIREMENTS)
	$(PIP_PATH) install --upgrade black flake8 pylint isort pytest pytest-benchmark

test:
	@echo "$(BLUE)Running tests...$(RESET)"
	$(PYTHON_PATH) -m pytest $(TEST_DIR) -v -s -x -m "not benchmark"

bench:
	@echo "$(BLUE)Running benchmarks...$(RESET)"
	$(PYTHON_PATH) -m pytest $(BENCH_DIR) -v -s -x -m "benchmark" --benchmark-columns=min,max,mean,stddev,median,rounds,outliers --benchmark-sort=name --benchmark-min-rounds=10

lint:
	@echo "$(BLUE)Running linters...$(RESET)"
	$(PYTHON_PATH) -m flake8 $(SRC_DIR) $(TEST_DIR)
	$(PYTHON_PATH) -m pylint $(SRC_DIR) $(TEST_DIR)

format:
	@echo "$(BLUE)Formatting code...$(RESET)"
	$(PYTHON_PATH) -m black $(SRC_DIR) $(TEST_DIR)
	$(PYTHON_PATH) -m isort $(SRC_DIR) $(TEST_DIR)

clean:
	@echo "$(BLUE)Cleaning project...$(RESET)"
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type d -name ".pytest_cache" -exec rm -r {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	rm -rf .coverage htmlcov
