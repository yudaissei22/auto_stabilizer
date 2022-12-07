#ifndef Stabilizer_H
#define Stabilizer_H

#include "GaitParam.h"
#include <prioritized_qp_osqp/prioritized_qp_osqp.h>
#include <cnoid/JointPath>

class Stabilizer{
public:
  // Stabilizerでしか使わないパラメータ
  std::vector<double> bodyAttitudeControlGain=std::vector<double>{0.5, 0.5}; // 要素数2. gaitParam.footMidCoords座標系X軸, Y軸. 単位は[/s]. 0以上
  std::vector<double> bodyAttitudeControlTimeConst=std::vector<double>{1000, 1000}; // 要素数2. gaitParam.footMidCoords座標系X軸, Y軸. 単位は[s]. 0より大きい
  std::vector<double> bodyAttitudeControlCompensationLimit=std::vector<double>{0.7,0.7}; //要素数2. gaitParam.footMidCoords座標系X軸, Y軸. 単位は[rad]. STが動いている間は変更されない. 0以上

  std::vector<std::vector<double> > supportPgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > supportDgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > landingPgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > landingDgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > swingPgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > swingDgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  double swing2LandingTransitionTime = 0.05; // [s]. 0より大きい
  double landing2SupportTransitionTime = 0.1; // [s]. 0より大きい
  double support2SwingTransitionTime = 0.2; // [s]. 0より大きい

  std::vector<cnoid::Vector6> ee_K; // 要素数EndEffectors. EndEffector frame. endEffector origin. 0以上
  std::vector<cnoid::Vector6> ee_D; // 要素数EndEffectors. EndEffector frame. endEffector origin. 0以上
  cnoid::Vector3 com_K; // generate frame. 
  cnoid::Vector3 com_D; // generate frame.
  cnoid::Vector6 root_K; // generate frame.

  void init(const GaitParam& gaitParam, cnoid::BodyPtr& actRobotTqc){
    for(int i=0;i<NUM_LEGS;i++){
      cnoid::JointPath jointPath(actRobotTqc->rootLink(), actRobotTqc->link(gaitParam.eeParentLink[i]));
      if(jointPath.numJoints() == 6){
        // supportPgain[i] = {5,10,10,5,0.1,0.1};
        // supportDgain[i] = {10,20,20,10,10,10};
        // landingPgain[i] = {5,10,10,5,0.1,0.1};
        // landingDgain[i] = {10,20,20,10,10,10};
        // swingPgain[i] = {5,10,10,5,0.1,0.1};
        // swingDgain[i] = {10,20,20,10,10,10};
        // 下はもとのauto_stabilizerの値. ゲインが低すぎて、go-velocity 0 0 0のときに前に進んでいってしまう
        /*supportPgain[i] = {5,15,10,5,0.2,0.2};
        supportDgain[i] = {10,20,20,10,5,5};
        landingPgain[i] = {5,15,1,1,0.2,0.2};
        landingDgain[i] = {10,10,10,10,5,5};
        swingPgain[i] = {5,30,20,10,5,5};
        swingDgain[i] = {10,30,20,20,30,30};*/
        supportPgain[i] = {0,0,0,0,0,0};
        supportDgain[i] = {0,0,0,0,0,0};
        landingPgain[i] = {0,0,0,0,0,0};
        landingDgain[i] = {0,0,0,0,0,0};
        swingPgain[i] = {0,0,0,0,0,0};
        swingDgain[i] = {0,0,0,0,0,0};
      }else{
        supportPgain[i].resize(jointPath.numJoints(), 100.0);
        supportDgain[i].resize(jointPath.numJoints(), 100.0);
        landingPgain[i].resize(jointPath.numJoints(), 100.0);
        landingDgain[i].resize(jointPath.numJoints(), 100.0);
        swingPgain[i].resize(jointPath.numJoints(), 100.0);
        swingDgain[i].resize(jointPath.numJoints(), 100.0);
      }
    }
    for(int i=0;i<gaitParam.eeName.size();i++){
      cnoid::Vector6 defaultEED;
      if(i<NUM_LEGS){
	defaultEED << 10, 10, 10, 10, 10, 10;
      }else{
	defaultEED << 10, 10, 10, 10, 10, 10;
      }
      this->ee_D.push_back(defaultEED);
      cnoid::Vector6 defaultEEK;
      if(i<NUM_LEGS){
	defaultEEK << 50, 50, 50, 20, 20, 20;
      }else{
	defaultEEK << 50, 50, 50, 20, 20, 20;
      }
      this->ee_K.push_back(defaultEEK);
    }
    cnoid::Vector3 defaultComK;
    defaultComK << 0, 0, 0;
    this->com_K = defaultComK;
    cnoid::Vector3 defaultComD;
    defaultComD << 0, 0, 0;
    this->com_D = defaultComD;
    cnoid::Vector6 defaultRootK;
    defaultRootK << 1, 1, 1, 1, 1, 1;
    this->root_K = defaultRootK;
  }
protected:
  // 計算高速化のためのキャッシュ. クリアしなくても別に副作用はない.
  // for calcWrench
  mutable std::shared_ptr<prioritized_qp_osqp::Task> constraintTask_ = std::make_shared<prioritized_qp_osqp::Task>();
  mutable std::shared_ptr<prioritized_qp_osqp::Task> tgtZmpTask_ = std::make_shared<prioritized_qp_osqp::Task>();
  mutable std::shared_ptr<prioritized_qp_osqp::Task> copTask_ = std::make_shared<prioritized_qp_osqp::Task>();
  // for calcTorque
  mutable std::shared_ptr<prioritized_qp_osqp::Task> eeComTask_ = std::make_shared<prioritized_qp_osqp::Task>();
  mutable std::shared_ptr<prioritized_qp_osqp::Task> rootTask_ = std::make_shared<prioritized_qp_osqp::Task>();
public:
  void initStabilizerOutput(const GaitParam& gaitParam,
                            cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Vector3& o_stTargetZmp, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage) const;

  bool execStabilizer(const GaitParam& gaitParam, double dt, bool useActState,
                      cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Position& o_stTargetRootPose, cnoid::Vector3& o_stTargetZmp, std::vector<cnoid::Vector6>& o_stEETargetWrench, cnoid::Vector3& o_stTargetCogAcc) const;

  bool calcTorque(double dt, const GaitParam& gaitParam, bool useActState, cnoid::BodyPtr& actRobotTqc, const std::vector<cnoid::Vector6>& tgtEEWrench /* 要素数EndEffector数. generate座標系. EndEffector origin*/,
                  std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage, Eigen::VectorXd& prev_q, Eigen::VectorXd& prev_dq, std::vector<cnoid::Vector6>& eePoseDiff_prev, std::vector<cnoid::Position>& eeTargetPosed, std::vector<cnoid::Position>& eeTargetPosedd) const;
protected:
  bool moveBasePosRotForBodyRPYControl(double dt, const GaitParam& gaitParam, bool useActState,
                                       cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Position& o_stTargetRootPose) const;
  bool calcZMP(const GaitParam& gaitParam, double dt, bool useActState,
               cnoid::Vector3& o_tgtZmp/*generate座標系*/, cnoid::Vector3& o_tgtForce/*generate座標系*/, cnoid::Vector3& o_tgtCogAcc /*generate座標系*/) const;
  bool calcWrench(const GaitParam& gaitParam, const cnoid::Vector3& tgtZmp/*generate座標系*/, const cnoid::Vector3& tgtForce/*generate座標系. ロボットが受ける力*/, bool useActState,
                  std::vector<cnoid::Vector6>& o_tgtEEWrench /* 要素数EndEffector数. generate座標系. EndEffector origin*/) const;
};

#endif
