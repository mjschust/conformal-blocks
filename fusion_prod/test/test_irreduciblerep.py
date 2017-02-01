'''
Created on Nov 19, 2016

@author: mjschust
'''
import unittest
import fusion_prod.cbd as cbd

class Test(unittest.TestCase):


    def testSL2IrrRep(self):
        
        liealg = cbd.TypeALieAlgebra(1)
        wt = cbd._Weight(liealg, [0])
        rep = cbd.IrrRep(liealg, wt)
        self.assertEqual(1, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = cbd._Weight(liealg, [1])
        rep = cbd.IrrRep(liealg, wt)
        self.assertEqual(2, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = cbd._Weight(liealg, [2])
        rep = cbd.IrrRep(liealg, wt)
        self.assertEqual(3, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [2])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        char_wt = cbd._Weight(liealg, [0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
    def testSL3IrrRep(self):
        liealg = cbd.TypeALieAlgebra(2)
        rep = cbd.IrrRep(liealg, [0,0])
        self.assertEqual(1, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")

        rep = cbd.IrrRep(liealg, [1,0])
        self.assertEqual(3, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [1, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")

        rep = cbd.IrrRep(liealg, [0,1])
        self.assertEqual(3, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [1, 0])
        self.assertFalse(char_wt in dom_char, "Character incorrect")
        char_wt = cbd._Weight(liealg, [0, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")

        rep = cbd.IrrRep(liealg, [1,1])
        self.assertEqual(8, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [1, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        char_wt = cbd._Weight(liealg, [0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(2, dom_char[char_wt], "Character incorrect")

        rep = cbd.IrrRep(liealg, [2,1])
        self.assertEqual(15, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [2, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        char_wt = cbd._Weight(liealg, [1, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(2, dom_char[char_wt], "Character incorrect")
        char_wt = cbd._Weight(liealg, [0, 2])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
    def testSL4IrrRep(self):
        liealg = cbd.TypeALieAlgebra(3)
        wt = cbd._Weight(liealg, [0, 0, 0])
        rep = cbd.IrrRep(liealg, wt)
        self.assertEqual(1, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [0, 0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = cbd._Weight(liealg, [1, 0, 0])
        rep = cbd.IrrRep(liealg, wt)
        self.assertEqual(4, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [1, 0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = cbd._Weight(liealg, [0, 0, 1])
        rep = cbd.IrrRep(liealg, wt)
        self.assertEqual(4, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [1, 0, 0])
        self.assertFalse(char_wt in dom_char, "Character incorrect")
        char_wt = cbd._Weight(liealg, [0, 0, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = cbd._Weight(liealg, [0, 1, 0])
        rep = cbd.IrrRep(liealg, wt)
        self.assertEqual(6, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [0, 1, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        
        wt = cbd._Weight(liealg, [1, 1, 1])
        rep = cbd.IrrRep(liealg, wt)
        self.assertEqual(64, rep.get_dimension(), "Dimension not correct")
        dom_char = rep.get_dominant_character()
        char_wt = cbd._Weight(liealg, [1, 1, 1])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(1, dom_char[char_wt], "Character incorrect")
        char_wt = cbd._Weight(liealg, [2, 0, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(2, dom_char[char_wt], "Character incorrect")
        char_wt = cbd._Weight(liealg, [0, 1, 0])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(4, dom_char[char_wt], "Character incorrect")
        char_wt = cbd._Weight(liealg, [0, 0, 2])
        self.assertTrue(char_wt in dom_char, "Character incorrect")
        self.assertEqual(2, dom_char[char_wt], "Character incorrect")
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSL2IrrRep', 'testSL3IrrRep']
    unittest.main()