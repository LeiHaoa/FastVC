#ifndef DATA_POOL_H
#define DATA_POOL_H

#include <vector>
#include "../Variation.h"
class dataPool{
public:
	vector<Variation*> _data;
	int _pointer;
	int _size;

	dataPool(int init_size){
		//_data.resize(init_size, new Variation());
		_data.reserve(init_size);
		_size = init_size;
		for(int i = 0; i < init_size; i++){
			_data.emplace_back(new Variation());
		}
	}
	void reset(){
		_pointer = 0;
		for(Variation* v: _data){
			memset(v, 0, sizeof(Variation));
		}
	}
	Variation* get_variation(){
		if(_pointer == (_size - 1)){
			//_data.emplace_back(Variation());
			_data.reserve(_size + _size/2);
			for(int i = _size; i < _data.capacity(); i++){
				_data.emplace_back(new Variation());
			}
			_size = _size + _size/2;
			//cerr << "resize data buf to: " << _size << "capacity: " << _data.capacity() << endl;
		}
		Variation* var = _data[_pointer];
		_pointer++;
		return var;
	}
};

#endif
