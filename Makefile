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
SRC_DIR = gte
REQUIREMENTS = requirements.txt
DEV_REQUIREMENTS = requirements-dev.txt

# Colors
GREEN = \033[0;32m
NC = \033[0m

# Targets
.PHONY: help venv venv-activate install install-dev test lint format clean run

help:
	@echo "Available commands:"
	@echo "  make venv           - Create virtual environment"
	@echo "  make activate       - Activate virtual environment (prints command)"
	@echo "  make install        - Install production dependencies"
	@echo "  make install-dev    - Install development dependencies"
	@echo "  make test           - Run tests"
	@echo "  make lint           - Run linters (flake8, pylint)"
	@echo "  make format         - Format code (black, isort)"
	@echo "  make clean          - Clean project"
	@echo "  make run            - Run main script"

venv:
	@echo "$(GREEN)Creating virtual environment...$(NC)"
	$(PYTHON) -m venv $(VENV_DIR)
	@echo "$(GREEN)Virtual environment created in $(VENV_DIR)$(NC)"
	@echo "To activate, run:"
	@echo "  source $(VENV_ACTIVATE)"

activate:
	@echo "$(GREEN)Run this command to activate virtual environment:$(NC)"
	@echo "source $(VENV_ACTIVATE)"

install:
	@echo "$(GREEN)Installing production dependencies...$(NC)"
	$(PIP_PATH) install --upgrade -r $(REQUIREMENTS)

install-dev: install
	@echo "$(GREEN)Installing development dependencies...$(NC)"
	$(PIP_PATH) install --upgrade -r $(DEV_REQUIREMENTS)
	$(PIP_PATH) install --upgrade black flake8 pylint isort pytest

test:
	@echo "$(GREEN)Running tests...$(NC)"
	$(PYTHON_PATH) -m pytest $(TEST_DIR) -v -s

lint:
	@echo "$(GREEN)Running linters...$(NC)"
	$(PYTHON_PATH) -m flake8 $(SRC_DIR) $(TEST_DIR)
	$(PYTHON_PATH) -m pylint $(SRC_DIR) $(TEST_DIR)

format:
	@echo "$(GREEN)Formatting code...$(NC)"
	$(PYTHON_PATH) -m black $(SRC_DIR) $(TEST_DIR)
	$(PYTHON_PATH) -m isort $(SRC_DIR) $(TEST_DIR)

clean:
	@echo "$(GREEN)Cleaning project...$(NC)"
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type d -name ".pytest_cache" -exec rm -r {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	rm -rf .coverage htmlcov

run:
	@echo "$(GREEN)Running project...$(NC)"
	$(PYTHON_PATH) $(SRC_DIR)/main.py