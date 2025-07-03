# HeatVizBackend/urls.py

from django.contrib import admin
from django.urls import include, path, re_path
from drf_yasg import openapi
from drf_yasg.views import get_schema_view
from rest_framework import permissions
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.reverse import reverse
from django.urls import get_resolver

# Schema view configuration
schema_view = get_schema_view(
    openapi.Info(
        title="HeatViz API",
        default_version='v1',
        description="API for handling heatmap data",
        contact=openapi.Contact(email="your_email@example.com"),
        license=openapi.License(name="BSD License"),
    ),
    public=True,
    permission_classes=(permissions.AllowAny,),
)

# API root view
@api_view(['GET'])
def api_root(request, format=None):
    return Response({
        'heatmapdata': reverse('process-data', request=request, format=format),
        'cleanup-session': reverse('cleanup-session', request=request, format=format),
        'command-execution': reverse('command-execution', request=request, format=format),
    })

urlpatterns = [
    path('', api_root),
    path('admin/', admin.site.urls),
    path('api/heatmapdata/', include('heatviz.urls')),
    # Schema paths served with Swagger/Redoc
    re_path(r'^swagger(?P<format>\.json|\.yaml)$', schema_view.without_ui(cache_timeout=0), name='schema-json'),
    re_path(r'^swagger/$', schema_view.with_ui('swagger', cache_timeout=0), name='schema-swagger-ui'),
    re_path(r'^redoc/$', schema_view.with_ui('redoc', cache_timeout=0), name='schema-redoc'),
]


# Print registered URL patterns for debugging
print("Registered URL patterns:")
for url_pattern in get_resolver().url_patterns:
    print(url_pattern)