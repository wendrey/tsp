/*
 * SampleDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

/***************************
 * Modified: 01.12.2016	   *
 * Wendrey Lustosa Cardoso *
 * RA: 148234			   *
 ***************************/

#include "SampleDecoder.h"

SampleDecoder::SampleDecoder(const LpdTspInstance &inst) {
	l = &inst;
}

SampleDecoder::~SampleDecoder() { }

// Dado um cromossomo, devolve o custo do circuito gerado pelo mesmo
// Atribui penalidades se encontrar alguma inviabilidade no circuito

double SampleDecoder::decode(const std::vector< double >& chromosome) const {
	
	int k = 0;
	double cost = 0.0;
	int penalty = 10000;
	
	// Dado um cromossomo, gera uma sequência de vértices
	
	std::vector <std::pair<double,DNode>> ranking(chromosome.size());

	cout << "cromossomo: " ;
	for (DNodeIt n(l->g); n != INVALID; ++n) {
		if (n != l->depot) { 
			ranking[k] = std::pair <double,DNode> (chromosome[k],n);
			k++;
			cout << chromosome[k];
		}
	}
	
	vector <DNode> tour;
	tour.push_back(l->depot);

	cout << " cromossomo: ";

	for(std::vector<std::pair<double,DNode>>::const_iterator i = ranking.begin(); i != ranking.end(); ++i)	
{		tour.push_back(i->second);
		cout << l->vname[i->second];	} 
		cout << endl;

	// Verifica se existe aresta ligando os vertices

	for(int i = 0; i < (int) tour.size(); i++) {
		OutArcIt o(l->g, tour[i]);
		for(; o != INVALID; ++o)
			if(l->g.target(o) == tour[(i+1) % (int)tour.size()])
				break;
		if(o == INVALID)
			cost += penalty;
		else
			cost += l->weight[o];
	}
	
	// Verifica a ordem de coleta e entrega dos itens
	
	vector<bool> ss(l->k, false);
	vector<bool> tt(l->k, false);

	for(int i = 1; i < (int)tour.size(); i++) {
		if(l->s[tour[i]] > 0) {
			if(!ss[l->s[tour[i]]-1] && !tt[l->s[tour[i]]-1])
				ss[l->s[tour[i]]-1] = true;
			else
				cost += penalty;
		}
		else if(l->t[tour[i]] > 0) {
			if(ss[l->t[tour[i]]-1] && !tt[l->t[tour[i]]-1])
				tt[l->t[tour[i]]-1] = true;
			else
				cost += penalty;
		}
	}

	// Verifica se a quantidade carregada é menor que a capacidade permitida

	double load = 0.0;

	for(int i = 0; i < (int)tour.size(); i++) {
		if(l->t[tour[i]] > 0 )
			load -= l->items[l->t[tour[i]]-1].w;
		if(l->s[tour[i]] > 0 )
			load += l->items[l->s[tour[i]]-1].w;
		if(load < (-1)*MY_EPS)
			cost += penalty;
		if(load > l->capacity)
			cost += penalty;
	}

	// Retorna o custo do circuito
	
	return cost;

}

