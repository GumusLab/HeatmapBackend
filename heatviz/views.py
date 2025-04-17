


# import os
# from tempfile import NamedTemporaryFile
# from rest_framework.decorators import api_view, parser_classes
# from rest_framework.parsers import MultiPartParser
# from rest_framework.response import Response
# from rest_framework import status
# from django.http import HttpResponse
# from heatviz.clustergrammer2.make_clustergrammer import make_cluster


# @api_view(['POST'])
# @parser_classes([MultiPartParser])
# def process_data_view(request):
#     file = request.FILES.get('data')  # Get the file from the request
#     print("*** file coming here is as follows ********", file)
#     if not file:
#         return Response({"error": "No file provided"}, status=status.HTTP_400_BAD_REQUEST)

#     try:
#         # Save the file to a temporary location
#         with NamedTemporaryFile(delete=False, suffix=".tsv") as temp_file:
#             for chunk in file.chunks():
#                 temp_file.write(chunk)
#             temp_file_path = temp_file.name

#         # Call the make_cluster function with the file path
#         response_data = make_cluster(temp_file_path)

#         # Create HTTP response
#         response = HttpResponse(response_data, content_type="application/json", status=200)
#         return response

#     except Exception as e:
#         return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

#     finally:
#         # Clean up the temporary file
#         if os.path.exists(temp_file_path):
#             os.remove(temp_file_path)


import os
import traceback
import logging
from tempfile import NamedTemporaryFile
from rest_framework.decorators import api_view, parser_classes
from rest_framework.parsers import MultiPartParser
from rest_framework.response import Response
from rest_framework import status
from django.http import HttpResponse, JsonResponse
from heatviz.clustergrammer2.make_clustergrammer import make_cluster

# Set up logger
logger = logging.getLogger('heatviz')

@api_view(['POST'])
@parser_classes([MultiPartParser])
def process_data_view(request):
    file = request.FILES.get('data')  # Get the file from the request
    print("*** file coming here is as follows ********", file)
    if not file:
        return Response({"error": "No file provided"}, status=status.HTTP_400_BAD_REQUEST)

    temp_file_path = None
    try:
        # Save the file to a temporary location
        with NamedTemporaryFile(delete=False, suffix=".tsv") as temp_file:
            for chunk in file.chunks():
                temp_file.write(chunk)
            temp_file_path = temp_file.name
        
        print(f"File saved to temporary location: {temp_file_path}")
        
        # Call the make_cluster function with the file path
        response_data = make_cluster(temp_file_path)

        # Write the JSON data to the file
        with open("top_250.json", 'w') as json_file:
            json_file.write(response_data)
        
        
        # Create HTTP response
        response = HttpResponse(response_data, content_type="application/json", status=200)
        return response

    except Exception as e:
        # Get detailed error information
        error_trace = traceback.format_exc()
        error_msg = f"Error processing file: {str(e)}"
        
        # Log the full error
        print(f"ERROR: {error_msg}")
        print(f"TRACEBACK: {error_trace}")
        
        # Log to file if logger is configured
        logger.error(f"{error_msg}\n{error_trace}")
        
        # Return error response with details
        return JsonResponse({
            "error": error_msg,
            "traceback": error_trace
        }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    finally:
        # Clean up the temporary file
        if temp_file_path and os.path.exists(temp_file_path):
            try:
                os.remove(temp_file_path)
                print(f"Temporary file {temp_file_path} removed")
            except Exception as e:
                print(f"Error removing temporary file: {e}")
