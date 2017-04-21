'''
Created on Nov 19, 2016

@author: mjschust
'''
import unittest
import conformal_blocks.cbbundle as cbd

class Test(unittest.TestCase):


    def testSL2IrrRep(self):
        
        liealg = cbd.TypeALieAlgebra(1)
        wt = tuple([0])
        self.assertEqual(1, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = tuple([1])
        self.assertEqual(2, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = tuple([2])
        self.assertEqual(3, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([2])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        char_wt = tuple([0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
    def testSL3IrrRep(self):
        liealg = cbd.TypeALieAlgebra(2)
        wt= (0,0)
        self.assertEqual(1, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")

        wt = (1, 0)
        self.assertEqual(3, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([1, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")

        wt = (0, 1)
        self.assertEqual(3, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([1, 0])
        self.assertFalse(char_wt in dom_char, "Character incorrect")
        char_wt = tuple([0, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")

        wt = (1, 1)
        self.assertEqual(8, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([1, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        char_wt = tuple([0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(2, dom_char[char_wt], "Character incorrect")

        wt = (2, 1)
        self.assertEqual(15, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([2, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        char_wt = tuple([1, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(2, dom_char[char_wt], "Character incorrect")
        char_wt = tuple([0, 2])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
    def testSL4IrrRep(self):
        liealg = cbd.TypeALieAlgebra(3)
        wt = tuple([0, 0, 0])
        self.assertEqual(1, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([0, 0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = tuple([1, 0, 0])
        self.assertEqual(4, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([1, 0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = tuple([0, 0, 1])
        self.assertEqual(4, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([1, 0, 0])
        self.assertFalse(char_wt in dom_char, "Character incorrect")
        char_wt = tuple([0, 0, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = tuple([0, 1, 0])
        self.assertEqual(6, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([0, 1, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = tuple([1, 1, 1])
        self.assertEqual(64, liealg.get_rep_dim(wt), "Dimension not correct")
        dom_char = liealg.get_dominant_character(wt)
        char_wt = tuple([1, 1, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        char_wt = tuple([2, 0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(2, dom_char[char_wt], "Character incorrect")
        char_wt = tuple([0, 1, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(4, dom_char[char_wt], "Character incorrect")
        char_wt = tuple([0, 0, 2])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(2, dom_char[char_wt], "Character incorrect")
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSL2IrrRep', 'testSL3IrrRep']
    unittest.main()