/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "StringUtils.h"
#include "Ciphertext.h"
#include "Scheme.h"

//----------------------------------------------------------------------------------
//   SHOW ARRAY
//----------------------------------------------------------------------------------

void StringUtils::show(uint64_t* vals, long size) {
	cout << "[";
	for (long i = 0; i < size; ++i) {
		cout << vals[i] << ", ";
	}
	cout << "]" << endl;
}

void StringUtils::show(long* vals, long size) {
	cout << "[";
	for (long i = 0; i < size; ++i) {
		cout << vals[i] << ", ";
	}
	cout << "]" << endl;
}

void StringUtils::show(double* vals, long size) {
	cout << "[";
	for (long i = 0; i < size; ++i) {
		cout << vals[i] << ", ";
	}
	cout << "]" << endl;
}

void StringUtils::show(complex<double>* vals, long size) {
	cout << "[";
	for (long i = 0; i < size; ++i) {
		cout << vals[i] << ", ";
	}
	cout << "]" << endl;
}


//----------------------------------------------------------------------------------
//   SHOW & COMPARE ARRAY
//----------------------------------------------------------------------------------


void StringUtils::showcompare(double val1, double val2, string prefix) {
	cout << "---------------------" << endl;
	cout << "m" + prefix + ":" << val1 << endl;
	cout << "d" + prefix + ":" << val2 << endl;
	cout << "e" + prefix + ":" << val1-val2 << endl;
	cout << "---------------------" << endl;
}

void StringUtils::showcompare(complex<double> val1, complex<double> val2, string prefix) {
	cout << "---------------------" << endl;
	cout << "m" + prefix + ":" << val1 << endl;
	cout << "d" + prefix + ":" << val2 << endl;
	cout << "e" + prefix + ":" << val1-val2 << endl;
	cout << "---------------------" << endl;
}

void StringUtils::showcompare(double* vals1, double* vals2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;
		cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;
		cout << "e" + prefix + ": " << i << " :" << (vals1[i]-vals2[i]) << endl;
		cout << "---------------------" << endl;
	}
}

void StringUtils::showcompare(complex<double>* vals1, complex<double>* vals2, long size, string prefix, bool PRT) {
    if(PRT==1){
        int correct=1;
        //Original
        for (long i = 0; i < size; ++i) {
            cout << "---------------------" << endl;
            cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;
            cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;
            cout << "e" + prefix + ": " << i << " :" << (vals1[i]-vals2[i]) << endl;
            cout << "---------------------" << endl;
            if( vals1[i].real()-vals2[i].real()>1e-5||
                vals1[i].real()-vals2[i].real()<-1e-5||
                vals1[i].imag()-vals2[i].imag()>1e-5||
                vals1[i].imag()-vals2[i].imag()<-1e-5)
                correct=0;
        }
        if (correct){
            cout<<prefix+"\t CORRECT!"<<endl<<endl;
        }
    }
    else{
        int correct=1;
        for (long i = 0; i < size; ++i) {
            if(vals1[i].real()-vals2[i].real()>1e-5||
                    vals1[i].real()-vals2[i].real()<-1e-5
            ){
                cout << "---------------------" << endl;
                cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;
                cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;
                cout << "e" + prefix + ": " << i << " :" << (vals1[i]-vals2[i]) << endl;
                cout << "---------------------" << endl;
                correct=0;
            }
        }
        if (correct){
            cout<<prefix+"\t CORRECT!"<<endl<<endl;
        }
    }
}


void StringUtils::showcompare(double* vals1, double val2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;
		cout << "d" + prefix + ": " << i << " :" << val2 << endl;
		cout << "e" + prefix + ": " << i << " :" << vals1[i]-val2 << endl;
		cout << "---------------------" << endl;
	}
}

void StringUtils::showcompare(complex<double>* vals1, complex<double> val2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;
		cout << "d" + prefix + ": " << i << " :" << val2 << endl;
		cout << "e" + prefix + ": " << i << " :" << vals1[i]-val2 << endl;
		cout << "---------------------" << endl;
	}
}

void StringUtils::showcompare(double val1, double* vals2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << val1 << endl;
		cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;
		cout << "e" + prefix + ": " << i << " :" << val1-vals2[i] << endl;
		cout << "---------------------" << endl;
	}
}

void StringUtils::showcompare(complex<double> val1, complex<double>* vals2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << val1 << endl;
		cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;
		cout << "e" + prefix + ": " << i << " :" << val1-vals2[i] << endl;
		cout << "---------------------" << endl;
	}
}

void StringUtils::showcompare(uint64_t* vals1, uint64_t* vals2, long size, string prefix) {
    for (long i = 0; i < size; ++i) {
        cout << "---------------------" << endl;
        cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;
        cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;
        cout << "e" + prefix + ": " << i << " :" << (vals1[i]-vals2[i]) << endl;
        cout << "---------------------" << endl;
    }
}

void StringUtils::CompareAndSet(long LimbsNum, long N, uint64_t *val, uint64_t *golden_val, string STRING, bool ComparePRT, bool Set)
{
    cout<< "Compare " + STRING <<endl;
    int AllCorrect=1;
    for (long i = 0; i < LimbsNum; ++i) {
        uint64_t *axj=val+(i*N);
        uint64_t *gvj=golden_val+(i*N);
        for (int j = 0; j < N; ++j) {
            //Compare After_ModRaise [Eval]
            if (axj[j]!=gvj[j]){
                if (ComparePRT)
                    cout<<i<<"\t"<<j<<endl;
                AllCorrect=0;
            }
            //SET After_ModRaise_Eval
            if(Set)
                axj[j]=gvj[j];
        }
    }
    if(AllCorrect==1) {
        cout<< STRING + " is the same! " <<endl;
    }
    else{
        cout<< STRING + " DOES NOT MATCH THE GOLDEN VALUE "<<endl;
        if (Set)
            cout<<"[!!] input is set to golden value "<<endl;
    }
    cout<<endl;
}

void StringUtils::CompareAndSet(long LimbsNum, long N, Ciphertext& ct, uint64_t *golden_val, string STRING, bool ComparePRT, bool Set)
{
    cout<< "Compare " + STRING <<endl;
    int AllCorrectAX=1;
    for (long i = 0; i < LimbsNum; ++i) {
        uint64_t *axj=ct.ax+(i*N);
        uint64_t *gvj=golden_val+(i+LimbsNum)*N;
        for (int j = 0; j < N; ++j) {
            //Compare After_ModRaise [Eval]
            if (axj[j]!=gvj[j]){
                if (ComparePRT)
                    cout<<"ax: "<<i<<"\t"<<j<<endl;
                AllCorrectAX=0;
            }
            //SET After_ModRaise_Eval
            if(Set)
                axj[j]=gvj[j];
        }
    }
    if(AllCorrectAX==1) {
        cout<< STRING + " Ax is the same! " <<endl;
    }
    else{
        cout<< STRING + " Ax DOES NOT MATCH THE GOLDEN VALUE "<<endl;
        if (Set)
            cout<<"[!!] Ax input is set to golden value "<<endl;
    }
    int AllCorrectBX=1;
    for (long i = 0; i < LimbsNum; ++i) {
        uint64_t *bxj=ct.bx+(i*N);
        uint64_t *gvbj=golden_val+(i*N);
        for (int j = 0; j < N; ++j) {
            //Compare After_ModRaise [Eval]
            if (bxj[j]!=gvbj[j]){
                if (ComparePRT)
                    cout<<i<<"\t"<<j<<endl;
                AllCorrectBX=0;
            }
            //SET After_ModRaise_Eval
            if(Set)
                bxj[j]=gvbj[j];
        }
    }
    if(AllCorrectBX==1) {
        cout<< STRING + " Bx is the same! " <<endl;
    }
    else{
        cout<< STRING + " Bx DOES NOT MATCH THE GOLDEN VALUE "<<endl;
        if (Set)
            cout<<"[!!] Bx input is set to golden value "<<endl;
    }
    cout<<endl;
}

void StringUtils::DecAndShow(Scheme scheme, Ciphertext ct, string String)
{
    cout<<"Decrypt "+String+"\t";
    complex<double> *dvecModRaiseRes = scheme.decrypt(scheme.secretKey, ct);
    StringUtils::show(dvecModRaiseRes, ct.slots); cout<<endl;
}

void StringUtils::FullDecAndShow(Scheme scheme, Ciphertext ct, string String)
{
    cout<<"Decrypt "+String+"\t";
    complex<double> *dvecModRaiseRes = scheme.fulldecrypt(scheme.secretKey, ct);
    StringUtils::show(dvecModRaiseRes, ct.slots); cout<<endl;
}
