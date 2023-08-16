import unittest
import importlib
from tempfile import TemporaryDirectory
import io
from contextlib import redirect_stdout

from merge_vcfs import Genotype, Variant, insert_variant, genotype_from_string, combine_sorted_lists, combine_variants

mendelian_consistency = importlib.import_module('mendelian-consistency')

class TestGenotype(unittest.TestCase):

	def test_phased(self):
		genotype = Genotype([0,1,2], True)
		self.assertEqual(genotype.get_alleles(), [0,1,2])
		self.assertEqual(str(genotype), "0|1|2")
		self.assertEqual(genotype.is_phased(), True)
		self.assertEqual(genotype.get_ploidy(), 3)

	def test_unphased(self):
		genotype = Genotype([0,1,2,3], False)
		self.assertEqual(genotype.get_alleles(), [0,1,2,3])
		self.assertEqual(str(genotype), "0/1/2/3")
		self.assertEqual(genotype.is_phased(), False)
		self.assertEqual(genotype.get_ploidy(), 4)


class TestVariant(unittest.TestCase):

	def test_init(self):
		samples = ['sample1', 'sample2']
		start = 100
		end = 101
		ref_allele = 'A'
		alt_alleles = ['T', 'AT']
		genotypes = [Genotype([0,1], True), Genotype([1,2], True)]
		ploidy = 2

		variant = Variant(samples, start, ref_allele, alt_alleles, genotypes, ploidy)
		self.assertEqual(variant.get_samples(), samples)
		self.assertEqual(variant.get_start(), start)
		self.assertEqual(variant.get_end(), end)
		self.assertEqual(variant.get_ref(), ref_allele)
		self.assertEqual(variant.get_alt(), alt_alleles)
		self.assertEqual(variant.get_ploidy(), ploidy)

		self.assertEqual(str(variant.get_genotype_of('sample1')), '0|1')
		self.assertEqual(str(variant.get_genotype_of('sample2')), '1|2')

		self.assertEqual(variant.get_allele(0), ref_allele)
		self.assertEqual(variant.get_allele(1), alt_alleles[0])
		self.assertEqual(variant.get_allele(2),alt_alleles[1])

	def test_operators(self):
		v1 = Variant(['sample1'], 200, 'AT', ['A'], [Genotype([1,1], True)], 2)
		v2 = Variant(['sample2'], 201, 'T', ['A'], [Genotype([1,1], True)], 2)
		v3 = Variant(['sample2'], 200, 'T', ['A'], [Genotype([1,1], True)], 2)

		self.assertEqual(v1 < v2, True)
		self.assertEqual(v1 <= v2, True)
		self.assertEqual(v1 > v2, False)
		self.assertEqual(v1 >= v2, False)

		self.assertEqual(v1 < v3, False)
		self.assertEqual(v1 <= v3, True)
		self.assertEqual(v1 > v3, False)
		self.assertEqual(v1 >= v3, True)


class TestInsertVariant(unittest.TestCase):

	def test_insert1(self):
		v = Variant(['sample1'], 16, 'C', ['CAAA'], [Genotype([1,1], True)], 2)
		self.assertEqual(v.get_allele(1), 'CAAA')
		sequence = 'ATCCGATGGCATTAGC'
		result, offset = insert_variant(sequence, v, 1, 13)
		self.assertEqual(result, 'ATCCAAAGATGGCATTAGC')
		self.assertEqual(offset, 10)

	def test_insert2(self):
		v = Variant(['sample1'], 16, 'C', ['CAAA'], [Genotype([1,1], True)], 2)
		self.assertEqual(v.get_allele(1), 'CAAA')
		sequence = 'ATCCGATGGCATTAGC'
		result, offset = insert_variant(sequence, v, 1, 13)
		self.assertEqual(result, 'ATCCAAAGATGGCATTAGC')
		self.assertEqual(offset, 10)

		v2 = Variant(['sample1'], 21, 'G', ['GCCC'], [Genotype([1,1], True)], 2)
		self.assertEqual(v2.get_allele(1), 'GCCC')
		result2, offset2 = insert_variant(result, v2, 1, offset)
		self.assertEqual(result2, 'ATCCAAAGATGGCCCCATTAGC')
		self.assertEqual(offset2, 7)

	def test_insert3(self):
		v1 = Variant(['sample1'], 24, 'A', ['ACCCC'], [Genotype([1,0],True)], 2)
		v2 = Variant(['sample1'], 33, 'AAA', ['T'], [Genotype([1,1], True)], 2)
		v3 = Variant(['sample1'], 38, 'A', ['C'], [Genotype([0,0], True)], 2)
		sequence = 'AAAAAAAAAAAAAAAAA'
		offset = 23
		sequence, offset = insert_variant(sequence, v1, 1, offset)
		self.assertEqual(sequence, 'AACCCCAAAAAAAAAAAAAAA')
		self.assertEqual(offset, 19)

		sequence, offset = insert_variant(sequence, v2, 1, offset)
		self.assertEqual(sequence, 'AACCCCAAAAAAAATAAAA')
		self.assertEqual(offset, 21)

		sequence, offset = insert_variant(sequence, v3, 1, offset)
		self.assertEqual(sequence, 'AACCCCAAAAAAAATAACA')
		self.assertEqual(offset, 21)

	def test_insert4(self):
		v = Variant(['sample1'], 0, 'AAAA', ['TTTT'], [Genotype([1,0],True)], 2)
		sequence = 'AAAAAAAAA'
		sequence, offset = insert_variant(sequence, v, 1, 0)
		self.assertEqual(sequence, 'TTTTAAAAA')
		self.assertEqual(offset, 0)


class TestGenotypeFromString(unittest.TestCase):

	def test_init(self):
		genotype = genotype_from_string('1|0')
		self.assertEqual(str(genotype), '1|0')
		genotype = genotype_from_string('1|1')
		self.assertEqual(str(genotype), '1|1')
		genotype = genotype_from_string('0/1')
		self.assertEqual(str(genotype), '0/1')


class TestCombineSortedLists(unittest.TestCase):

	def test_combine(self):
		l1 = [1,2,3,4,5,6]
		l2 = [4,5,6,7]
		self.assertEqual(combine_sorted_lists(l1,l2), [1,2,3,4,4,5,5,6,6,7])

		l1 = [1,2,5,7]
		l2 = [0,4,8]
		self.assertEqual(combine_sorted_lists(l1,l2), [0,1,2,4,5,7,8])

class TestCombineVariants(unittest.TestCase):

	def test_simple(self):
		v1 = Variant(['sample1'], 2, 'AAA', ['TTT'], [Genotype([0,1], True)], 2)
		v2 = Variant(['sample2'], 7, 'AAA', ['CCC'], [Genotype([1,0], True)], 2)
		combined = combine_variants([v1,v2], {'chr1':'AAAAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 2, 'AAAAAAAA', ['TTTAAAAA', 'AAAAACCC'], [Genotype([0,1], True), Genotype([2,0], True)], 2)
		print('simple', combined, expected)
		self.assertEqual(combined, expected)

	def test_simple2(self):
		v1 = Variant(['sample1'], 2, 'AAA', ['TTT'], [Genotype([1,1], True)], 2)
		v2 = Variant(['sample2'], 7, 'AAA', ['CCC'], [Genotype([1,0], True)], 2)
		combined = combine_variants([v1,v2], {'chr1':'AAAAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 2, 'AAAAAAAA', ['TTTAAAAA', 'AAAAACCC'], [Genotype([1,1], True), Genotype([2,0], True)], 2)
		print('simple2', combined, expected)
		self.assertEqual(combined, expected)

	def test_simple3(self):
		v1 = Variant(['sample1'], 2, 'A', ['TTT'], [Genotype([1,1], True)], 2)
		v2 = Variant(['sample1', 'sample2'], 7, 'AAA', ['AC', 'AT'], [Genotype([0,1], True), Genotype([1,0], True)], 2)
		combined = combine_variants([v1,v2], {'chr1':'AAAAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 2, 'AAAAAAAA', ['TTTAAAAAAA', 'TTTAAAAAC', 'AAAAAAC'], [Genotype([1,2], True), Genotype([3,0], True)], 2)
		print('simple3', combined, expected)
		self.assertEqual(combined, expected)

	def test_nested1(self):
		v1 = Variant(['sample2'], 1, 'AAAAAAA', ['CCCCCCC'], [Genotype([1,0], True)], 2)
		v2 = Variant(['sample1'], 6, 'AAAAAAAAAAA', ['TTTTTTTTTTT'], [Genotype([1,0], True)], 2)
		v3 = Variant(['sample2'], 14, 'A', ['G'], [Genotype([1,1], True)], 2)
		combined = combine_variants([v1,v2, v3], {'chr1': 'AAAAAAAAAAAAAAAAAAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 1, 'AAAAAAAAAAAAAAAA', ['AAAAATTTTTTTTTTT', 'CCCCCCCAAAAAAGAA', 'AAAAAAAAAAAAAGAA', ], [Genotype([1,0], True), Genotype([2,3], True)], 2)
		print('nested1', combined, expected)
		self.assertEqual(combined, expected)

	def test_nested2(self):
		v1 = Variant(['sample1', 'sample2'], 1, 'AAA', ['CCC'], [Genotype([1,0,0], True), Genotype([0,1,0], True)], 3)
		v2 = Variant(['sample1'], 6, 'AAA', ['TTT'], [Genotype([0,1,0], True)], 3)
		v3 = Variant(['sample2'], 3, 'AAAA', ['GGGG'], [Genotype([1,0,0], True)], 3)
		v4 = Variant(['sample2'], 8, 'AA', ['GG'], [Genotype([0,1,0], True)], 3)
		combined = combine_variants([v1,v2,v3,v4], {'chr1':'AAAAAAAAAAAAAAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 1, 'AAAAAAAAA', ['CCCAAAAAA', 'AAAAATTTA', 'AAGGGGAAA', 'CCCAAAAGG'], [Genotype([1,2,0], True), Genotype([3,4,0], True)], 3)

		print('nested2', combined, expected)
		self.assertEqual(combined, expected)

	def test_nested3(self):
		v1 = Variant(['sample1'], 3, 'AAAA', ['TTTT'], [Genotype([1], True)], 1)
		v2 = Variant(['sample2'], 3, 'AA', ['GG'], [Genotype([1], True)], 1)
		v3 = Variant(['sample2'], 5, 'AA', ['CC'], [Genotype([1], True)], 1)
		combined = combine_variants([v1,v2,v3], {'chr1': 'AAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 3, 'AAAA', ['TTTT', 'GGCC'], [Genotype([1], True), Genotype([2], True)], 1)
		print('nested3', combined, expected)
		self.assertEqual(combined, expected)

	def test_overlapping(self):
		v1 = Variant(['sample1'], 3, 'AAAA', ['TTTT'], [Genotype([1], True)], 1)
		v2 = Variant(['sample1'], 4, 'AAAA', ['CCCC'], [Genotype([1], True)], 1)

		combined = combine_variants([v1,v2], {'chr1': 'AAAAAAAAAA'}, 'chr1')
		self.assertEqual(combined, None)

		v3 = Variant(['sample1'], 4, 'AAAA', ['CCCC'], [Genotype([0], True)], 1)
		combined = combine_variants([v1,v3], {'chr1': 'AAAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1'], 3, 'AAAA', ['TTTT'], [Genotype([1], True)], 1)

	def test_one_absent(self):
		v1 = Variant(['sample1'], 3, 'AAAA', ['GGGG'], [Genotype([0,0], True)], 2)
		v2 = Variant(['sample2'], 5, 'AAAA', ['CCCC'], [Genotype([0,1], True)], 2)
		combined = combine_variants([v1,v2], {'chr1': 'AAAAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 5, 'AAAA', ['CCCC'], [Genotype([0,0], True), Genotype([0,1], True)], 2)
		print('one_absent', combined)
		self.assertEqual(combined, expected)

	def test_all_absent(self):
		v1 = Variant(['sample1'], 3, 'AAAA', ['GGGG'], [Genotype([0,0], True)], 2)
		v2 = Variant(['sample2'], 5, 'AAAA', ['CCCC'], [Genotype([0,0], True)], 2)
		combined = combine_variants([v1,v2], {'chr1': 'AAAAAAAAAAA'}, 'chr1')
		expected = None
		self.assertEqual(combined, expected)

	def test_undefined(self):
		v1 = Variant(['sample1'], 3, 'AAAA', ['TTTT'], [Genotype([1], True)], 1)
		v2 = Variant(['sample2'], 3, 'AA', ['GG'], [Genotype([1], True)], 1)
		v3 = Variant(['sample2'], 5, 'AA', ['CC'], [Genotype([-1], True)], 1)
		combined = combine_variants([v1,v2,v3], {'chr1': 'AAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 3, 'AAAA', ['TTTT'], [Genotype([1], True), Genotype([-1], True)], 1)
		print('nested3', combined, expected)
		self.assertEqual(combined, expected)

	def test_undefined2(self):
		v1 = Variant(['sample2'], 1, 'AAAAAAA', ['CCCCCCC'], [Genotype([1,0], True)], 2)
		v2 = Variant(['sample1'], 6, 'AAAAAAAAAAA', ['TTTTTTTTTTT'], [Genotype([1,-1], True)], 2)
		v3 = Variant(['sample2'], 14, 'A', ['G'], [Genotype([1,1], True)], 2)
		combined = combine_variants([v1,v2, v3], {'chr1': 'AAAAAAAAAAAAAAAAAAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 1, 'AAAAAAAAAAAAAAAA', ['AAAAATTTTTTTTTTT', 'CCCCCCCAAAAAAGAA', 'AAAAAAAAAAAAAGAA', ], [Genotype([1,-1], True), Genotype([2,3], True)], 2)
		print('nested1', combined, expected)
		self.assertEqual(combined, expected)

	def test_undefined3(self):
		v1 = Variant(['sample1'], 3, 'AAAA', ['TTTT'], [Genotype([-1], True)], 1)
		v2 = Variant(['sample2'], 3, 'AA', ['GG'], [Genotype([1], True)], 1)
		v3 = Variant(['sample2'], 5, 'AA', ['CC'], [Genotype([1], True)], 1)
		combined = combine_variants([v1,v2,v3], {'chr1': 'AAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 3, 'AAAA', ['GGCC'], [Genotype([-1], True), Genotype([1], True)], 1)
		print('nested3', combined, expected)
		self.assertEqual(combined, expected)

	def test_undefined4(self):
		v1 = Variant(['sample1'], 3, 'AAAA', ['TTTT'], [Genotype([-1], True)], 1)
		v2 = Variant(['sample2'], 3, 'AA', ['GG'], [Genotype([-1], True)], 1)
		v3 = Variant(['sample2'], 5, 'AA', ['CC'], [Genotype([-1], True)], 1)
		combined = combine_variants([v1,v2,v3], {'chr1': 'AAAAAAAAA'}, 'chr1')
		expected = Variant(['sample1', 'sample2'], 3, 'AAAA', ['GGCC'], [Genotype([-1], True), Genotype([1], True)], 1)
		print('nested3', combined, expected)
		self.assertEqual(combined, expected)


#class TestCheckMendelianConsistency(unittest.TestCase):
#
#	def test_cons1(self):
#		child = Genotype([0,1], True)
#		mother = Genotype([1,1], False)
#		father = Genotype([1,0], True)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, mother, father), True)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, father, mother), True)
#		child = Genotype([0,0], False)
#		mother = Genotype([0,1], True)
#		father = Genotype([0,1], False)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, mother, father), True)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, father, mother), True)

#	def test_incons1(self):
#		child = Genotype([1,1], True)
#		mother = Genotype([0,0], False)
#		father = Genotype([1,1], False)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, mother, father), False)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, father, mother), False)
#		child = Genotype([0,2], True)
#		mother = Genotype([1,1], True)
#		father = Genotype([0,2], True)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, mother, father), False)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, father, mother), False)
#		child = Genotype([0,4], False)
#		mother = Genotype([1,4], False)
#		father = Genotype([3,4], False)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, mother, father), False)
#		self.assertEqual(mendelian_consistency.check_mendelian_consistency(child, father, mother), False)


#class TestRemoveSamplesFromHeader(unittest.TestCase):
#
#	def test_remove_samples(self):
#		header='CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\tsample4'
#		samples=['sample1', 'sample3']
#		updated = mendelian_consistency.remove_samples_from_header(header, samples)
#		expected='CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample2\tsample4'
#		self.assertEqual(updated, expected)
#		samples=[]
#		updated = mendelian_consistency.remove_samples_from_header(header, samples)
#		self.assertEqual(updated, header)
#		samples = ['sample4']
#		updated = mendelian_consistency.remove_samples_from_header(header, samples)
#		expected = 'CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3'
#		self.assertEqual(updated, expected)

#class TestRunMendelian(unittest.TestCase):
#	
#	def test_example(self):
#		vcffile='testfiles/test.vcf'
#		pedfile='testfiles/test.ped'

#		with TemporaryDirectory() as tempdir:
#			for (file, r) in [('testfiles/expected-all.vcf', False), ('testfiles/expected-parents.vcf', True)]:
#				outstats = tempdir + '/outstats.csv'
#				outvcf = tempdir + '/outvcf.vcf'

#				f = io.StringIO()
#				with redirect_stdout(f):
#					mendelian_consistency.run_mendelian(vcffile, pedfile, outstats, r)
#				outfile = open(outvcf, 'w')
#				outfile.write(f.getvalue())
#				outfile.close()
#				print(f.getvalue())

				# compare to expected vcf
#				line_count = 0
#				for computed, expected in zip(open(outvcf, 'r'), open(file, 'r')):
#					line_count += 1
#					computed_fields = computed.split()
#					expected_fields = expected.split()
#					self.assertEqual(computed_fields, expected_fields)
#				self.assertEqual(line_count, 9)
				# check computed statistics
#				for line in open(outstats, 'r'):
#					if line.startswith("#"):
#						continue
#					fields = line.split()
#					if fields[0] == 'HG00733':
#						self.assertEqual(fields[1], '9.0')
#						self.assertEqual(fields[2], '12.0')
#						continue
#					if fields[0] == 'NA19240':
#						self.assertEqual(fields[1], '8.0')
#						self.assertEqual(fields[2], '12.0')
#						continue
#					assert False

#	def test_example2(self):
#		vcffile='testfiles/test2.vcf'
#		pedfile='testfiles/test2.ped'

#		with TemporaryDirectory() as tempdir:
#			for (file, r) in [('testfiles/expected-all2.vcf', False), ('testfiles/expected-parents2.vcf', True)]:
#				outstats = tempdir + '/outstats.csv'
#				outvcf = tempdir + '/outvcf.vcf'

#				f = io.StringIO()
#				with redirect_stdout(f):
#					mendelian_consistency.run_mendelian(vcffile, pedfile, outstats, r)
#				outfile = open(outvcf, 'w')
#				outfile.write(f.getvalue())
#				outfile.close()
#				print(f.getvalue())

				# compare to expected vcf
#				line_count = 0
#				for computed, expected in zip(open(outvcf, 'r'), open(file, 'r')):
#					line_count += 1
#					computed_fields = computed.split()
#					expected_fields = expected.split()
#					self.assertEqual(computed_fields, expected_fields)
#				self.assertEqual(line_count, 8)
				# check computed statistics
#				for line in open(outstats, 'r'):
#					if line.startswith("#"):
#						continue
#					fields = line.split()
#					if fields[0] == 'HG00733':
#						self.assertEqual(fields[1], '9.0')
#						self.assertEqual(fields[2], '12.0')
#						continue
#					if fields[0] == 'NA19240':
#						self.assertEqual(fields[1], '8.0')
#						self.assertEqual(fields[2], '12.0')
#						continue
#					if fields[0] == 'HG00734':
#						self.assertEqual(fields[1], '10.0')
#						self.assertEqual(fields[2], '12.0')
#						continue
#					assert False
			

if __name__ == '__main__':
    unittest.main()


