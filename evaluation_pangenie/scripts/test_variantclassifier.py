import unittest
import importlib
from tempfile import TemporaryDirectory
import io
from contextlib import redirect_stdout
import vcf

from variantclassifier import determine_type_from_ids, VariantType, determine_variant_type, determine_pyvcf_type, determine_variant_length, determine_variant_class, determine_class_from_line

class Test_determine_type_from_ids(unittest.TestCase):

	def test_determine_type_from_ids1(self):
		ids = ['chr1-1-SNV-0-1', 'chr1-1-SNV-1-1', 'chr1-1-INS-1-50', 'chr1-1-DEL-1-49', 'chr1-1-COMPLEX-1-50', 'chr1-1-COMPLEX-1-1', 'chr1-1-COMPLEX-1-20']
		vartypes = determine_type_from_ids(ids)
		print(vartypes)
		expected = [VariantType.snp, VariantType.snp, VariantType.large_insertion, VariantType.midsize_deletion, VariantType.large_complex, VariantType.small_complex, VariantType.midsize_complex]
		
		self.assertEqual(set(vartypes), set(expected))
		
	def test_determine_type_from_ids2(self):
		var_id = ['chr1-1-SNV-0-1']
		vartype = determine_type_from_ids(var_id)
		assert len(vartype) == 1
		self.assertEqual(len(vartype), 1)
		self.assertEqual(vartype[0], VariantType.snp)

		var_id = ['chr1-1-INS-0-1']
		vartype = determine_type_from_ids(var_id)
		assert len(vartype) == 1
		self.assertEqual(len(vartype), 1)
		self.assertEqual(vartype[0], VariantType.small_insertion)
		
		var_id = ['chr1-1-DEL-0-19']
		vartype = determine_type_from_ids(var_id)
		assert len(vartype) == 1
		self.assertEqual(len(vartype), 1)
		self.assertEqual(vartype[0], VariantType.small_deletion)
		
		var_id = ['chr1-1-INS-0-20']
		vartype = determine_type_from_ids(var_id)
		assert len(vartype) == 1
		self.assertEqual(len(vartype), 1)
		self.assertEqual(vartype[0], VariantType.midsize_insertion)

		var_id = ['chr1-1-INS-0-50']
		vartype = determine_type_from_ids(var_id)
		assert len(vartype) == 1
		self.assertEqual(len(vartype), 1)
		self.assertEqual(vartype[0], VariantType.large_insertion)
		
class Test_determine_variant_type(unittest.TestCase):

	def test_determine_variant_type(self):
		expected = [
			VariantType.snp,
			VariantType.large_insertion,
			VariantType.snp,
			VariantType.large_deletion,
			VariantType.snp,
			VariantType.small_complex,
			VariantType.small_complex,
			VariantType.large_complex,
			VariantType.large_complex
		]
		
		computed = []
		for record in vcf.Reader(filename='data/small1.vcf'):
			computed.append(determine_variant_type(record))
		self.assertEqual(computed, expected)

class Test_determine_pyvcf_type(unittest.TestCase):

	def test_determine_pyvcf_type(self):
		expected = [
			VariantType.snp,
			VariantType.large_insertion,
			VariantType.snp,
			VariantType.large_deletion,
			VariantType.snp,
			VariantType.small_complex,
			VariantType.small_complex,
			VariantType.large_complex,
			VariantType.midsize_insertion
		]
		
		computed = []
		for record in vcf.Reader(filename='data/small1.vcf'):
			computed.append(determine_pyvcf_type(record))
		self.assertEqual(computed, expected)
		
class Test_determine_variant_length(unittest.TestCase):
	def test_determine_variant_length(self):
		expected = [1,50,1,50,1,2,2,50,50]
		computed = []
		for record in vcf.Reader(filename='data/small1.vcf'):
			computed.append(determine_variant_length(record))
		self.assertEqual(computed, expected)
		
class Test_determine_variant_class(unittest.TestCase):
	def test_determine_variant_class(self):
		expected = ['SNV', 'INS', 'SNV', 'DEL', 'SNV', 'COMPLEX', 'COMPLEX', 'COMPLEX', 'COMPLEX']
		computed = []
		for record in vcf.Reader(filename='data/small1.vcf'):
			computed.append(determine_variant_class(record))
		self.assertEqual(computed, expected)

class Test_determine_variant_class_from_line(unittest.TestCase):
	def test_determine_variant_class(self):
		expected = ['SNV', 'INS', 'SNV', 'DEL', 'SNV', 'COMPLEX', 'COMPLEX', 'COMPLEX', 'INS']
		computed = []
		for record in open('data/small1.vcf', 'r'):
			if record.startswith('#'):
				continue
			computed.append(determine_class_from_line(record))
		self.assertEqual(computed, expected)
