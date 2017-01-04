#pragma once
/* T is a data structure, such as an OctNode, to provide a problem space.
 * AIPROCESS is a Artificial Intelligent class to make decision.
 * \example
 *
 *
 */
#include <vector>

namespace SU
{
	enum ACT_FLG {
		FLG_REMOVE,
		FLG_ADD,
		FLG_UNCHANGE
	};

	/* A rule engine for agent
	*
	*/
	template<typename AgentType>
	class suRuleEngine
	{
	public:
		float cutRatio;
		double maxStrain;
		int boundaryThickness;
		suRuleEngine() {}
		ACT_FLG processStrain(AgentType *pAgent) {
			bool cut = false;
			if (pAgent->pState_->strainSort <=(1-pow(cutRatio,globalValue::globalValuePoint().iteratorTimes))&&pAgent->pState_->location>boundaryThickness&&pAgent->pState_->strain<=maxStrain)
				cut = true;
			if (cut == true)
				return FLG_REMOVE;
			else
				return FLG_UNCHANGE;
		}
		ACT_FLG process(AgentType *pAgent)
		{
			//std::cout << "suRuleEngine is running \n";
			bool cut = true;
			if (pAgent->pState_->label_ == SU::INTERIOR_CELL) {
				//std::cout << std::endl;
				//std::cout << pAgent->pState_->location;
				for (int i = 0; i < pAgent->environment_.size(); i++) {
					if (pAgent->environment_[i]->out==true&&pAgent->pState_->location < pAgent->environment_[i]->location) {
						//std::cout << pAgent->environment_[i]->location << std::endl;
						cut = false;
						break;
					}
				}
			}
			else
			{
				cut = false;
			}
			if (cut == true)
				return FLG_REMOVE;
			else
				return FLG_UNCHANGE;
			/*if (pAgent->pState_->location == 3)
				return FLG_REMOVE;
			else
				return FLG_UNCHANGE;*/
		}
	};

	/* Editor for OctNodes, which represents agent action
	*
	*/
	template<typename AgentType>
	class suEditor
	{
	public:		
		void process(int id, AgentType *pAgent)
		{
			if (id == FLG_REMOVE) {
				pAgent->kill();
				return;
			}
			if (id == FLG_UNCHANGE) {
				return;
			}
			if (id == FLG_ADD) {
				pAgent->grow();
				return;
			}
			
			//std::cout << "Editor process\n";
		}

	};

	template <typename T>
	class suAgent
	{
	public:
		suAgent() : isAlive_(true), pState_(0) {}

		int make_decision(float ratio_,double strain_,int thickness_)
		{
			//make decision by rules
			suRuleEngine<suAgent<OctNode> > ai;
			ai.cutRatio = ratio_;
			ai.maxStrain = strain_;
			ai.boundaryThickness = thickness_;
			return ai.processStrain(this);
		} 
		void act(float ratio,double strain,int thickness)
		{
			int act_flag = make_decision(ratio,strain,thickness);
			suEditor<suAgent<OctNode> > editor;
			editor.process(act_flag, this);
		}
		void set_environment(std::vector<T *> env)
		{
			environment_.clear();
			environment_ = env;
		}

		void update();
		void kill() { isAlive_ = false; }
		void grow() { isAlive_ = true; }
		bool alive() { return isAlive_; }
		void bind(T *pObject, std::vector<T *> env)//pObject is current agent    env is the neighbors
		{
			pState_ = pObject;
			environment_ = env;
		}
		T *pState_;
		std::vector<T *> environment_;
	private:
		bool isAlive_;
		
		T tempState_;
		
		std::vector<T> tempEnvironment_;
	};

	///implementation
	template<typename T>
	void suAgent<T>::update()
	{
		if (isAlive_ == true) {
			pState_->out = true;
		}
		else {
			pState_->out = false;
		}
		//std::cout << "Update()" << std::endl;
	}


	

}


