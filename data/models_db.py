from sqlalchemy import Column, Integer, String, DateTime
from datetime import datetime
from core.database import Base

class ApprovedVariant(Base):
    __tablename__ = "approved_variants"

    id = Column(Integer, primary_key=True, index=True)
    chromosome = Column(String, index=True)
    position = Column(Integer, index=True)
    ref_allele = Column(String)
    alt_allele = Column(String)
    gene = Column(String, nullable=True)
    approved_by = Column(String)
    created_at = Column(DateTime, default=datetime.utcnow)
    job_id = Column(String, nullable=True)
    patient_id = Column(String, nullable=True)
    assembly = Column(String, default="GRCh38")
