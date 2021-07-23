import unittest as ut
import datamol as dm

from medchem.filter import demerits


class Test_DemeritsFilter(ut.TestCase):

    test_config = {
        "output": "test",
        "min_atoms": 7,
        "soft_max_atoms": 30,
        "hard_max_atoms": 50,
        "smarts": [],
        "nodemerit": False,
        "dthresh": 160,
        "odm": [],
        "okiso": False,
        "noapdm": False,
    }
    smiles_list = [
        "Cc1cnc(CNc(cccc2-c3cn(CC(C4)CC4O)c4ncnc(N)c34)c2F)s1"
        "Cc1cnc(CNc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc(cccc2-c3cn(C[C@H](C4)C[C@H]4O)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc(cccc2-c3cn(C[C@H]4C[C@H](CO)C4)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc(cccc2-c3cn(C[C@H]4C[C@@H](CO)C4)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc2c(CO)c(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)ccc2)s1",
        "CC(C)(C)c1cnc(CNc(cccc2-c3cn(CC(C4)CC4O)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc2nc(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)ccc2)s1",
        "Cc1nc(CNc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)cs1",
        "NC(c(ccc(NC(c1cccc(C(F)(F)F)c1)=O)c1)c1-c1cc2cnc(NC3CC3)nc2cc1)=O",
        "Cc1cnc(CNc(cccc2-c3cn(Cc4cn(CCN5CCOCC5)nn4)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc2cc(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)ccc2)s1",
        "Cc1cnc(CNc(cccc2-c3cn(CC4CCOCC4)c4ncnc(N)c34)c2F)s1",
        "Cc1cnc(CNc(cccc2-c3cn(C[C@@H]4COCC4)c4ncnc(N)c34)c2F)s1",
        "Cn1nnc(CNc2c(c(-c3cccc(NCc4nccs4)c3F)c[nH]3)c3ncn2)c1",
        "Cn1c(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)cnc1",
        "Cc1cnc(C[C@@H]2NCc3c2cccc3-c2cn(C[C@H](C3)C[C@@H]3O)c3ncnc(N)c23)s1",
        "Cc1ncc(Cn2nnc(CNc3ncnc4c3c(-c(cccc3NCc5nccs5)c3F)c[nH]4)c2)s1",
        "CCc1noc(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)n1",
        "CCn1c(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)cnc1",
        "Fc1cccc(-n2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)c1",
        "CCc1nc(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)no1",
        "Cc(n1CCn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)ncc1[N+]([O-])=O",
        "Fc(cc1F)cc(Br)c1-n1nnc(CNc2ncnc3c2c(-c(cccc2NCc4nccs4)c2F)c[nH]3)c1",
        "Fc1cc(-n2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)cc(F)c1",
        "Cc1cnc(CNc(cc2)cc(-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)s1",
        "C[C@@H](c1cn(-c2cc(C(F)(F)F)cc(Br)c2)nn1)Nc1ncnc2c1c(-c(cccc1NCc3nccs3)c1F)c[nH]2",
        "Nc1c(c(-c2cccc(NCc3cscn3)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(Cc4cncc(F)c4)nn3)c2ncn1",
        "C[C@@H](c1cn(Cc2cncn2C2CC2)nn1)Nc1c(c(-c2cccc(NCc3nccs3)c2F)c[nH]2)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3cscn3)c2F)cn2C[C@H](C3)C[C@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(-c4cc([N+]([O-])=O)ccc4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C[C@@H]4OCCc5c4cccc5)nn3)c2ncn1",
        "CC(c1cn(CC2C(CC3)CC3C2)nn1)N(Cc1nccs1)c(cccc1-c2c[nH]c3ncnc(N)c23)c1F",
        "Nc1c(c(-c2cccc(NCc3ncc4n3CCC4)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "NC1CC(Cn2c3ncnc(N)c3c(-c3cccc(NCc4nc(CO)cs4)c3F)c2)C1",
        "Nc1c(c(-c2cccc(NCc3ncccc3)c2F)nn2CC(C3)CC3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CC5(COCC5)OCC4)nn3)c2ncn1",
        "C[C@@H](c1cn(CC2CCCCC2)nn1)N(Cc1nccs1)c(cccc1-c2c[nH]c3ncnc(N)c23)c1F",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CC4C(CC5)CC5C4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CC5(CSCC5)OCC4)nn3)c2ncn1",
        "CCn1ncnc1Cn1nnc(CNc2c(c(-c3cccc(NCc4nccs4)c3F)c[nH]3)c3ncn2)c1",
        "Nc1c(c(-c(cccc2NCc3nccs3)c2F)cn2Cc3cn([C@@H]4CC5(CCSCC5)OCC4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn([C@@H]4CC5(CCCCC5)OCC4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn([C@@H]4CC5(CCC5)OCC4)nn3)c2ncn1",
        "CCn1c(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)ncc1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn([C@@H]4CC5(CCOCC5)OCC4)nn3)c2ncn1",
        "CC1OC(Cn2nnc(CNc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)CC1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C[C@@H](C4)c5c4cccc5)nn3)c2ncn1",
        "C[C@@H](c1cn(Cc2nnnn2C2CC2)nn1)Nc1c(c(-c2cccc(NCc3nccs3)c2F)c[nH]2)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCN4CCOCC4)nn3)c2ncn1",
        "Cc1cnc(C(Nc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)=O)s1",
        "Nc1c(c(-c(cccc2NCc3nccs3)c2F)cn2Cc3cn(-c(cc4)cc5c4scc5)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCC4CC4)nn3)c2ncn1",
        "CC(C)Cn1ncnc1Cn1nnc(CNc2c(c(-c3cccc(NCc4nccs4)c3F)c[nH]3)c3ncn2)c1",
        "Nc1c(c(-c2cccc(NCc3nc(CO)cs3)c2F)cn2CC3CCC3)c2ncn1",
        "CC(c1cn(C2CC3(CCCC3)OCC2)nn1)Nc1c(c(-c2cccc(NCc3nccs3)c2F)c[nH]2)c2ncn1",
        "Nc1c(c(-c2cccc(NCC(C(NC=C3)=O)=C3c3ncccc3)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCOCC4)nn3)c2ncn1",
        "CCN1CC(Cn2nnc(C(C)Nc3c(c(-c4cccc(NCc5nccs5)c4F)c[nH]4)c4ncn3)c2)OCC1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(-c4cc(F)cc(F)c4)nn3)c2ncn1",
        "CC(c1cn(C2CC3(CCOCC3)OCC2)nn1)Nc1c(c(-c2cccc(NCc3nccs3)c2F)c[nH]2)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCCCCC4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCN4CCCC4)nn3)c2ncn1",
        "Nc1c(c(-c2cccc(NCC(C(NC=C3)=O)=C3c3ccccc3)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCCCC4)nn3)c2ncn1",
        "Cc(n1CCn2nnc(Cn3c4ncnc(N)c4c(-c4cccc(NCc5nccs5)c4F)c3)c2)ncc1[N+]([O-])=O",
        "Nc1c(c(-c2cccc(NCc3ncc(-c4ccccc4)s3)c2F)cn2C[C@H](C3)C[C@@H]3O)c2ncn1",
        "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCCC4)nn3)c2ncn1",
    ]

    def test_demerits(self):
        res = demerits.demerit_filter(self.smiles_list, **self.test_config)


if __name__ == "__main__":
    ut.main()
