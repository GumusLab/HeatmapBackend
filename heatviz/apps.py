# from django.apps import AppConfig


# class HeatvizConfig(AppConfig):
#     default_auto_field = 'django.db.models.BigAutoField'
#     name = 'heatviz'

from django.apps import AppConfig

class HeatvizConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'heatviz'
    
    def ready(self):
        # Initialize pathway database when Django starts
        try:
            from .services.pathway_system import initialize_pathway_database
            initialize_pathway_database("pathway_data")
            print("✅ Pathway database initialized successfully")
        except Exception as e:
            print(f"⚠️  Warning: Could not initialize pathway database: {e}")