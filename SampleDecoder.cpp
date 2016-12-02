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

// Runs in \Theta(n \log n):
double SampleDecoder::decode(const std::vector< double >& chromosome) const {

/*	std::vector< std::pair< double, unsigned > > ranking(chromosome.size());

	for(unsigned i = 0; i < chromosome.size(); ++i) {
		ranking[i] = std::pair< double, unsigned >(chromosome[i], i);
	}

	// Here we sort 'permutation', which will then produce a permutation of [n] in pair::second:
	std::sort(ranking.begin(), ranking.end());

	// permutation[i].second is in {0, ..., n - 1}; a permutation can be obtained as follows
	std::list< unsigned > permutation;
	for(std::vector< std::pair< double, unsigned > >::const_iterator i = ranking.begin();
			i != ranking.end(); ++i) {
		permutation.push_back(i->second);
	}

	// sample fitness is the first allele
	return chromosome.front();
*/

	// ----------
	
	int k = 0;
	double cost = 0.0;
	int penalty = 1000000;
	std::vector <std::pair<double,DNode>> ranking(chromosome.size());

	for (DNodeIt n(l->g); n != INVALID; ++n) {
		if (n != l->depot) { 
			ranking[k] = std::pair <double,DNode> (chromosome[k],n);
			k++;
		}
	}
	
	vector <DNode> tour;
	tour.push_back(l->depot);

	for(std::vector<std::pair<double,DNode>>::const_iterator i = ranking.begin(); i != ranking.end(); ++i)	
		tour.push_back(i->second);
	
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
			if(!ss[l->s[tour[i]]] && !tt[l->s[tour[i]]])
				ss[l->s[tour[i]]] = true;
			else
				cost += penalty;
		}
		else if(l->t[tour[i]] > 0) {
			if(ss[l->t[tour[i]]] && !tt[l->t[tour[i]]])
				tt[l->t[tour[i]]] = true;
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
	
	return cost;

}

