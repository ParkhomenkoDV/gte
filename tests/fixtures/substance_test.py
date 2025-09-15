import pytest
from substance import Substance


@pytest.fixture
def substance_factory():
    """Фабрика substance"""

    def create_user(name="Test"):
        return Substance(
            name=name,
        )

    return create_user
