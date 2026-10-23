"""Módulo de API FastAPI assíncrona para predição de resistência antimicrobiana e relatórios XAI (Seção 9)."""

from phyloresistance.api.app import app, criar_app, obter_registro_modelos

__all__ = ["app", "criar_app", "obter_registro_modelos"]

