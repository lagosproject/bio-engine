import os
import sys
import logging
from fastapi.testclient import TestClient

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import app
from services.reference import load_reference

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

client = TestClient(app)

def test_api_download(accession="NG_008866"):
    logger.info(f"Testing API-based download for {accession}...")
    
    # 1. Test the /check-reference endpoint
    logger.info("Calling /check-reference...")
    response = client.get(f"/check-reference?id={accession}")
    logger.info(f"Response: {response.status_code} - {response.json()}")
    
    # 2. Test the internal load_reference via the app context
    # This simulates what happens inside /create-job
    logger.info(f"Calling ref_service.load_reference({accession}) within app context...")
    try:
        from services import reference as ref_service
        path = ref_service.load_reference(accession)
        logger.info(f"Successfully loaded reference at: {path}")
        
        # Verify file exists and has content
        if os.path.exists(path):
            size = os.path.getsize(path)
            logger.info(f"File size: {size} bytes")
            if size > 0:
                print("\nAPI CONTEXT DOWNLOAD TEST PASSED")
            else:
                print("\nAPI CONTEXT DOWNLOAD TEST FAILED: File is empty")
        else:
            print(f"\nAPI CONTEXT DOWNLOAD TEST FAILED: File {path} not found")
            
    except Exception as e:
        logger.error(f"API context download failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    test_api_download()
